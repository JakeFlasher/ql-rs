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
mod domain_suite_30_tests {
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
mod domain_suite_32_tests {
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

mod evaluate_fabozzi_treasury_futures_ctd_and_switch_option_module {
    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[inline]
    fn smin<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b - (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[inline]
    fn smax<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[derive(Clone, Copy, Debug)]
    pub struct FabozziTreasuryFuturesCtdSwitchSpec<T: Real = f64> {
        pub bond_coupons: [T; 2],
        pub bond_maturities: [T; 2],
        pub bond_clean_prices: [T; 2],
        pub bond_modified_durations: [T; 2],
        pub futures_price: T,
        pub delivery_horizon_years: T,
        pub repo_rate: T,
        pub yield_volatility: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct FabozziTreasuryFuturesCtdSwitchReport<T: Real = f64> {
        pub conversion_factors: [T; 2],
        pub gross_basis: [T; 2],
        pub net_basis_bnoc: [T; 2],
        pub implied_repo_rates: [T; 2],
        pub ctd_implied_repo_rate: T,
        pub ctd_net_basis_bnoc: T,
        pub quality_switch_option_value: T,
        pub switch_adjusted_theoretical_futures: T,
    }

    pub fn evaluate_fabozzi_treasury_futures_ctd_and_switch_option<T: Real>(
        spec: &FabozziTreasuryFuturesCtdSwitchSpec<T>,
    ) -> FabozziTreasuryFuturesCtdSwitchReport<T> {
        let mut cfs = [T::zero(); 2];
        let mut gbasis = [T::zero(); 2];
        let mut nbasis = [T::zero(); 2];
        let mut irrs = [T::zero(); 2];
        let mut f_conv = [T::zero(); 2];
        let tau = spec.delivery_horizon_years.max_val(T::from_f64(1e-6));
        for k in 0..2 {
            let disc = T::from_f64(1.03).pow(-T::from_f64(2.0) * spec.bond_maturities[k]);
            let cf = (spec.bond_coupons[k] / T::from_f64(0.06)) * (T::one() - disc) + disc;
            let p = spec.bond_clean_prices[k];
            let gb = p - spec.futures_price * cf;
            let carry = (T::from_f64(100.0) * spec.bond_coupons[k] - p * spec.repo_rate) * tau;
            let nb = gb - carry;
            let irr = (spec.futures_price * cf + T::from_f64(100.0) * spec.bond_coupons[k] * tau
                - p)
                / (p * tau);
            cfs[k] = cf;
            gbasis[k] = gb;
            nbasis[k] = nb;
            irrs[k] = irr;
            f_conv[k] = (p - carry) / cf;
        }
        let ctd_irr = smax(irrs[0], irrs[1]);
        let ctd_bnoc = smin(nbasis[0], nbasis[1]);
        let dv0 = spec.bond_modified_durations[0] * spec.bond_clean_prices[0] / cfs[0];
        let dv1 = spec.bond_modified_durations[1] * spec.bond_clean_prices[1] / cfs[1];
        let d_dv = dv1 - dv0;
        let sigma_s = ((d_dv * d_dv + T::from_f64(1e-6))
            * spec.yield_volatility
            * spec.yield_volatility
            * tau)
            .sqrt();
        let diff_f = f_conv[1] - f_conv[0];
        let d = -((diff_f * diff_f + T::from_f64(1e-12)).sqrt()) / sigma_s;
        let switch_opt = (-spec.repo_rate * tau).exp()
            * (sigma_s * norm_pdf(d) + d * sigma_s * norm_cdf(d)).max_val(T::from_f64(1e-8));
        let min_conv = smin(f_conv[0], f_conv[1]);
        FabozziTreasuryFuturesCtdSwitchReport {
            conversion_factors: cfs,
            gross_basis: gbasis,
            net_basis_bnoc: nbasis,
            implied_repo_rates: irrs,
            ctd_implied_repo_rate: ctd_irr,
            ctd_net_basis_bnoc: ctd_bnoc,
            quality_switch_option_value: switch_opt,
            switch_adjusted_theoretical_futures: min_conv - switch_opt,
        }
    }

    #[derive(Clone, Copy, Debug)]
    pub struct AyacheForsythVetzalConvertibleSpec<T: Real = f64> {
        pub spot_price: T,
        pub bond_par: T,
        pub coupon_rate: T,
        pub conversion_ratio: T,
        pub maturity_years: T,
        pub risk_free_rate: T,
        pub dividend_yield: T,
        pub equity_volatility: T,
        pub hazard_rate: T,
        pub bond_recovery_rate: T,
        pub equity_default_jump_loss_eta: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct AyacheForsythVetzalConvertibleReport<T: Real = f64> {
        pub tf_cash_only_bond_floor: T,
        pub tf_convertible_fair_value: T,
        pub afv_jump_compensated_equity_drift: T,
        pub afv_expected_default_recovery_pv: T,
        pub afv_survival_continuation_pv: T,
        pub afv_convertible_fair_value: T,
        pub tf_vs_afv_splitting_wedge: T,
    }

    fn bs_call<T: Real>(s: T, k: T, sigma: T, r: T, q: T, t: T) -> T {
        let st = sigma * t.sqrt();
        let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * sigma * sigma) * t) / st;
        let d2 = d1 - st;
        s * (-q * t).exp() * norm_cdf(d1) - k * (-r * t).exp() * norm_cdf(d2)
    }

    pub fn evaluate_ayache_forsyth_vetzal_default_coupled_cb<T: Real>(
        spec: &AyacheForsythVetzalConvertibleSpec<T>,
    ) -> AyacheForsythVetzalConvertibleReport<T> {
        let sc = (T::one() - spec.bond_recovery_rate) * spec.hazard_rate;
        let rc = (spec.risk_free_rate + sc).max_val(T::from_f64(1e-8));
        let t = spec.maturity_years.max_val(T::from_f64(1e-6));
        let b_tf = spec.bond_par * (-rc * t).exp()
            + spec.coupon_rate * spec.bond_par * (T::one() - (-rc * t).exp()) / rc;
        let k_tf = b_tf / spec.conversion_ratio;
        let v_tf = b_tf
            + spec.conversion_ratio
                * bs_call(
                    spec.spot_price,
                    k_tf,
                    spec.equity_volatility,
                    spec.risk_free_rate,
                    spec.dividend_yield,
                    t,
                );
        let drift_afv = spec.risk_free_rate - spec.dividend_yield
            + spec.hazard_rate * spec.equity_default_jump_loss_eta;
        let q_afv = spec.dividend_yield - spec.hazard_rate * spec.equity_default_jump_loss_eta;
        let r_rf = spec.risk_free_rate.max_val(T::from_f64(1e-8));
        let b_rf = spec.bond_par * (-r_rf * t).exp()
            + spec.coupon_rate * spec.bond_par * (T::one() - (-r_rf * t).exp()) / r_rf;
        let surv = (-spec.hazard_rate * t).exp()
            * (b_rf
                + spec.conversion_ratio
                    * bs_call(
                        spec.spot_price,
                        b_rf / spec.conversion_ratio,
                        spec.equity_volatility,
                        spec.risk_free_rate,
                        q_afv,
                        t,
                    ));
        let rl = (spec.risk_free_rate + spec.hazard_rate).max_val(T::from_f64(1e-8));
        let def_claim = smax(
            spec.bond_recovery_rate * spec.bond_par,
            spec.conversion_ratio
                * (T::one() - spec.equity_default_jump_loss_eta)
                * spec.spot_price,
        );
        let rec_pv = (spec.hazard_rate / rl) * (T::one() - (-rl * t).exp()) * def_claim;
        let v_afv = surv + rec_pv;
        let wedge_diff = v_afv - v_tf;
        AyacheForsythVetzalConvertibleReport {
            tf_cash_only_bond_floor: b_tf,
            tf_convertible_fair_value: v_tf,
            afv_jump_compensated_equity_drift: drift_afv,
            afv_expected_default_recovery_pv: rec_pv,
            afv_survival_continuation_pv: surv,
            afv_convertible_fair_value: v_afv,
            tf_vs_afv_splitting_wedge: (wedge_diff * wedge_diff + T::from_f64(1e-24)).sqrt(),
        }
    }

    #[cfg(test)]
    mod domain_suite_47_tests {
        use super::*;
        use ql_ad::Dual;
        fn ctd_spec<T: Real>(x: T) -> FabozziTreasuryFuturesCtdSwitchSpec<T> {
            FabozziTreasuryFuturesCtdSwitchSpec {
                bond_coupons: [T::from_f64(0.06), T::from_f64(0.045)],
                bond_maturities: [T::from_f64(7.0), T::from_f64(9.5)],
                bond_clean_prices: [T::from_f64(101.0) + x, T::from_f64(94.5)],
                bond_modified_durations: [T::from_f64(5.8), T::from_f64(7.6)],
                futures_price: T::from_f64(100.2),
                delivery_horizon_years: T::from_f64(0.25),
                repo_rate: T::from_f64(0.04),
                yield_volatility: T::from_f64(0.01),
            }
        }
        fn afv_spec<T: Real>(x: T) -> AyacheForsythVetzalConvertibleSpec<T> {
            AyacheForsythVetzalConvertibleSpec {
                spot_price: T::from_f64(100.0) + x,
                bond_par: T::from_f64(100.0),
                coupon_rate: T::from_f64(0.04),
                conversion_ratio: T::one(),
                maturity_years: T::from_f64(3.0),
                risk_free_rate: T::from_f64(0.03),
                dividend_yield: T::from_f64(0.015),
                equity_volatility: T::from_f64(0.25),
                hazard_rate: T::from_f64(0.025),
                bond_recovery_rate: T::from_f64(0.4),
                equity_default_jump_loss_eta: T::from_f64(0.7),
            }
        }
        #[test]
        fn suite47_fixed_income_invariants_and_dual_fd() {
            let c0 = evaluate_fabozzi_treasury_futures_ctd_and_switch_option(&ctd_spec(0.0));
            assert!((c0.conversion_factors[0] - 1.0).abs() < 1e-12);
            assert!(c0.quality_switch_option_value > 0.0);
            let h = 1e-6;
            let fd_c = (evaluate_fabozzi_treasury_futures_ctd_and_switch_option(&ctd_spec(h))
                .switch_adjusted_theoretical_futures
                - evaluate_fabozzi_treasury_futures_ctd_and_switch_option(&ctd_spec(-h))
                    .switch_adjusted_theoretical_futures)
                / (2.0 * h);
            let ad_c = evaluate_fabozzi_treasury_futures_ctd_and_switch_option(&ctd_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .switch_adjusted_theoretical_futures;
            assert!(ad_c.d[0].abs() > 1e-6 && (ad_c.d[0] - fd_c).abs() < 1e-5);

            let a0 = evaluate_ayache_forsyth_vetzal_default_coupled_cb(&afv_spec(0.0));
            assert!(a0.afv_jump_compensated_equity_drift > 0.03 - 0.015);
            assert!(a0.afv_convertible_fair_value >= a0.tf_cash_only_bond_floor);
            let fd_a = (evaluate_ayache_forsyth_vetzal_default_coupled_cb(&afv_spec(h))
                .afv_convertible_fair_value
                - evaluate_ayache_forsyth_vetzal_default_coupled_cb(&afv_spec(-h))
                    .afv_convertible_fair_value)
                / (2.0 * h);
            let ad_a = evaluate_ayache_forsyth_vetzal_default_coupled_cb(&afv_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .afv_convertible_fair_value;
            assert!(ad_a.d[0].abs() > 1e-6 && (ad_a.d[0] - fd_a).abs() < 1e-5);
        }
    }
}
pub use evaluate_fabozzi_treasury_futures_ctd_and_switch_option_module::*;

mod evaluate_fabozzi_callable_putable_oas_and_convexity_module {
    #![allow(clippy::needless_range_loop, clippy::assign_op_pattern)]
    use ql_ad::Real;

    #[inline]
    fn smin<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b - (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[inline]
    fn smax<T: Real>(a: T, b: T) -> T {
        let d = a - b;
        T::from_f64(0.5) * (a + b + (d * d + T::from_f64(1e-16)).sqrt())
    }

    #[derive(Clone, Copy, Debug)]
    pub struct FabozziCallablePutableOasSpec<T: Real = f64> {
        pub initial_short_rate: T,
        pub rate_drift_per_year: T,
        pub short_rate_vol: T,
        pub oas_spread: T,
        pub coupon_rate: T,
        pub call_strike_price: T,
        pub put_strike_price: T,
        pub maturity_years: T,
        pub yield_shift_dy: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct FabozziCallablePutableOasReport<T: Real = f64> {
        pub straight_bond_fair_value: T,
        pub callable_bond_fair_value: T,
        pub putable_bond_fair_value: T,
        pub embedded_call_option_value: T,
        pub embedded_put_option_value: T,
        pub embedded_option_identity_error: T,
        pub straight_effective_duration: T,
        pub callable_effective_duration: T,
        pub straight_effective_convexity: T,
        pub callable_effective_convexity: T,
    }

    fn price_trees_shifted<T: Real>(
        spec: &FabozziCallablePutableOasSpec<T>,
        shift: T,
    ) -> (T, T, T) {
        let dt = spec.maturity_years.max_val(T::from_f64(1e-6)) / T::from_f64(3.0);
        let par = T::from_f64(100.0);
        let cpn = par * spec.coupon_rate * dt;
        let term_payoff = par + cpn;
        let mut vs = [term_payoff; 4];
        let mut vc = [term_payoff; 4];
        let mut vp = [term_payoff; 4];
        let half = T::from_f64(0.5);
        let sig_dt = spec.short_rate_vol * dt.sqrt();

        for step in (0..3).rev() {
            let t_f = T::from_f64(step as f64);
            for j in 0..=step {
                let two_j_minus_t = T::from_f64((2 * j as i32 - step as i32) as f64);
                let r_node = (spec.initial_short_rate + t_f * spec.rate_drift_per_year)
                    * (two_j_minus_t * sig_dt).exp()
                    + spec.oas_spread
                    + shift;
                let df = (-r_node * dt).exp();
                let cont_s = df * half * (vs[j + 1] + vs[j]);
                let cont_c = df * half * (vc[j + 1] + vc[j]);
                let cont_p = df * half * (vp[j + 1] + vp[j]);
                if step > 0 {
                    vs[j] = cont_s + cpn;
                    vc[j] = smin(cont_c, spec.call_strike_price) + cpn;
                    vp[j] = smax(cont_p, spec.put_strike_price) + cpn;
                } else {
                    vs[j] = cont_s;
                    vc[j] = cont_c;
                    vp[j] = cont_p;
                }
            }
        }
        (vs[0], vc[0], vp[0])
    }

    pub fn evaluate_fabozzi_callable_putable_oas_and_convexity<T: Real>(
        spec: &FabozziCallablePutableOasSpec<T>,
    ) -> FabozziCallablePutableOasReport<T> {
        let dy = spec.yield_shift_dy.max_val(T::from_f64(1e-5));
        let (vs_0, vc_0, vp_0) = price_trees_shifted(spec, T::zero());
        let (vs_up, vc_up, _) = price_trees_shifted(spec, dy);
        let (vs_dn, vc_dn, _) = price_trees_shifted(spec, -dy);

        let call_val = vs_0 - vc_0;
        let put_val = vp_0 - vs_0;
        let diff = (vc_0 + call_val) - (vp_0 - put_val);
        let id_err = (diff * diff + T::from_f64(1e-30)).sqrt();

        let two = T::from_f64(2.0);
        let d_eff_s = (vs_dn - vs_up) / (two * vs_0 * dy);
        let d_eff_c = (vc_dn - vc_up) / (two * vc_0 * dy);
        let c_eff_s = (vs_up + vs_dn - two * vs_0) / (vs_0 * dy * dy);
        let c_eff_c = (vc_up + vc_dn - two * vc_0) / (vc_0 * dy * dy);

        FabozziCallablePutableOasReport {
            straight_bond_fair_value: vs_0,
            callable_bond_fair_value: vc_0,
            putable_bond_fair_value: vp_0,
            embedded_call_option_value: call_val,
            embedded_put_option_value: put_val,
            embedded_option_identity_error: id_err,
            straight_effective_duration: d_eff_s,
            callable_effective_duration: d_eff_c,
            straight_effective_convexity: c_eff_s,
            callable_effective_convexity: c_eff_c,
        }
    }

    #[cfg(test)]
    mod domain_suite_48_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fabozzi_callable_putable_oas_invariants_and_ad_parity() {
            let spec = FabozziCallablePutableOasSpec {
                initial_short_rate: 0.035,
                rate_drift_per_year: 0.002,
                short_rate_vol: 0.18,
                oas_spread: 0.006,
                coupon_rate: 0.055,
                call_strike_price: 100.5,
                put_strike_price: 99.5,
                maturity_years: 3.0,
                yield_shift_dy: 0.0025,
            };
            let rep = evaluate_fabozzi_callable_putable_oas_and_convexity(&spec);
            assert!(rep.embedded_option_identity_error < 1e-12);
            assert!(rep.callable_effective_duration < rep.straight_effective_duration);
            assert!(rep.callable_effective_convexity < rep.straight_effective_convexity);

            let h = 1e-6;
            let r_ad = evaluate_fabozzi_callable_putable_oas_and_convexity(
                &FabozziCallablePutableOasSpec {
                    initial_short_rate: Dual::<1>::variable(spec.initial_short_rate, 0),
                    rate_drift_per_year: Dual::<1>::from_value(spec.rate_drift_per_year),
                    short_rate_vol: Dual::<1>::from_value(spec.short_rate_vol),
                    oas_spread: Dual::<1>::from_value(spec.oas_spread),
                    coupon_rate: Dual::<1>::from_value(spec.coupon_rate),
                    call_strike_price: Dual::<1>::from_value(spec.call_strike_price),
                    put_strike_price: Dual::<1>::from_value(spec.put_strike_price),
                    maturity_years: Dual::<1>::from_value(spec.maturity_years),
                    yield_shift_dy: Dual::<1>::from_value(spec.yield_shift_dy),
                },
            );
            let mut sp = spec;
            sp.initial_short_rate += h;
            let mut sm = spec;
            sm.initial_short_rate -= h;
            let _fd = (evaluate_fabozzi_callable_putable_oas_and_convexity(&sp)
                .callable_bond_fair_value
                - evaluate_fabozzi_callable_putable_oas_and_convexity(&sm)
                    .callable_bond_fair_value)
                / (2.0 * h);
            assert!((r_ad.callable_bond_fair_value.d[0] - _fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_fabozzi_callable_putable_oas_and_convexity_module::*;

mod evaluate_fabozzi_collared_frn_discount_margin_module {
    // Fabozzi (2021) & Tuckman (2011) Collared Floating-Rate Note (`FRN`) Discount Margin (`DM`),
    // Spread Duration vs. Rate Duration & Option-Adjusted Discount Margin (`OADM`) Engine.

    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    #[derive(Debug, Clone, Copy)]
    pub struct FabozziCollaredFrnDiscountMarginSpec<T: Real = f64> {
        pub par_face: T,
        pub index_rate_r: T,
        pub quoted_margin_qm: T,
        pub discount_margin_dm: T,
        pub cap_strike: T,
        pub floor_strike: T,
        pub rate_normal_vol_sigma: T,
        pub time_to_next_reset_tau1: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct FabozziCollaredFrnDiscountMarginReport<T: Real = f64> {
        pub straight_frn_price: T,
        pub par_reset_identity_error: T,
        pub embedded_cap_pv: T,
        pub embedded_floor_pv: T,
        pub collared_frn_price: T,
        pub spread_duration_dm: T,
        pub interest_rate_duration: T,
        pub option_adjusted_discount_margin_oadm: T,
    }

    fn bachelier_caplet<T: Real>(fwd: T, strike: T, std: T, is_cap: bool) -> T {
        let s = std.max_val(T::from_f64(1e-8));
        let diff = if is_cap { fwd - strike } else { strike - fwd };
        let d = diff / s;
        diff * norm_cdf(d) + s * norm_pdf(d)
    }

    fn straight_frn_pv<T: Real>(par: T, r: T, qm: T, dm: T, tau1: T) -> T {
        let n = 4usize;
        let mut pv = T::zero();
        let c_rate = r + qm;
        let disc_rate = (r + dm).max_val(T::from_f64(1e-4));
        for k in 1..=n {
            let tk = tau1 + T::from_f64((k - 1) as f64);
            let df = (T::one() + disc_rate).pow(-tk);
            let accrual = if k == 1 { tau1 } else { T::one() };
            let cf = if k == n {
                par * (T::one() + c_rate * accrual)
            } else {
                par * c_rate * accrual
            };
            pv += cf * df;
        }
        pv
    }

    pub fn evaluate_fabozzi_collared_frn_discount_margin<T: Real>(
        spec: &FabozziCollaredFrnDiscountMarginSpec<T>,
    ) -> FabozziCollaredFrnDiscountMarginReport<T> {
        let par = spec.par_face;
        let r = spec.index_rate_r;
        let qm = spec.quoted_margin_qm;
        let dm = spec.discount_margin_dm;
        let tau1 = spec.time_to_next_reset_tau1.max_val(T::from_f64(0.01));

        let v_straight = straight_frn_pv(par, r, qm, dm, tau1);
        // Exact Fabozzi/Tuckman FRN par-reset identity: at reset (tau1 = 1) with DM == QM, V == Par!
        let v_par_check = straight_frn_pv(par, r, qm, qm, T::one());
        let par_err = (v_par_check - par).abs();

        let fwd_coupon = r + qm;
        let disc_rate = (r + dm).max_val(T::from_f64(1e-4));
        let mut cap_pv = T::zero();
        let mut floor_pv = T::zero();
        let mut risky_dv01 = T::zero();

        for k in 1..=4usize {
            let tk = tau1 + T::from_f64((k - 1) as f64);
            let df = (T::one() + disc_rate).pow(-tk);
            let std = spec.rate_normal_vol_sigma * tk.sqrt();
            cap_pv += par * df * bachelier_caplet(fwd_coupon, spec.cap_strike, std, true);
            floor_pv += par * df * bachelier_caplet(fwd_coupon, spec.floor_strike, std, false);
            risky_dv01 += par * df;
        }

        // Collared FRN = Straight FRN + Long Floor - Short Cap:
        let v_collared = v_straight + floor_pv - cap_pv;

        // Analytical spread duration (w.r.t. DM) vs. interest-rate duration (w.r.t. parallel index shift):
        let bump = T::from_f64(1e-4);
        let sd_dm = (straight_frn_pv(par, r, qm, dm - bump, tau1)
            - straight_frn_pv(par, r, qm, dm + bump, tau1))
            / (T::from_f64(2.0) * bump * v_straight);
        let ir_dur = (straight_frn_pv(par, r - bump, qm, dm, tau1)
            - straight_frn_pv(par, r + bump, qm, dm, tau1))
            / (T::from_f64(2.0) * bump * v_straight);

        // Option-Adjusted Discount Margin (OADM):
        let oadm = dm + (floor_pv - cap_pv) / risky_dv01.max_val(T::from_f64(1e-6));

        FabozziCollaredFrnDiscountMarginReport {
            straight_frn_price: v_straight,
            par_reset_identity_error: par_err,
            embedded_cap_pv: cap_pv,
            embedded_floor_pv: floor_pv,
            collared_frn_price: v_collared,
            spread_duration_dm: sd_dm,
            interest_rate_duration: ir_dur.abs(),
            option_adjusted_discount_margin_oadm: oadm,
        }
    }

    #[cfg(test)]
    mod domain_suite_52_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fabozzi_collared_frn_discount_margin_invariants_and_ad_parity() {
            let spec = FabozziCollaredFrnDiscountMarginSpec {
                par_face: 100.0,
                index_rate_r: 0.042,
                quoted_margin_qm: 0.012,
                discount_margin_dm: 0.015,
                cap_strike: 0.065,
                floor_strike: 0.035,
                rate_normal_vol_sigma: 0.011,
                time_to_next_reset_tau1: 0.25,
            };
            let rep = evaluate_fabozzi_collared_frn_discount_margin(&spec);
            // Par-reset theorem & DM > QM => discount price < 100:
            assert!(rep.par_reset_identity_error < 1e-11);
            assert!(rep.straight_frn_price < 100.0);
            // Spread duration of an FRN is much larger than its interest-rate duration:
            assert!(rep.spread_duration_dm > rep.interest_rate_duration * 4.0);

            let h = 1e-6;
            let r_ad = evaluate_fabozzi_collared_frn_discount_margin(
                &FabozziCollaredFrnDiscountMarginSpec {
                    par_face: Dual::<1>::from_value(spec.par_face),
                    index_rate_r: Dual::<1>::from_value(spec.index_rate_r),
                    quoted_margin_qm: Dual::<1>::from_value(spec.quoted_margin_qm),
                    discount_margin_dm: Dual::<1>::variable(spec.discount_margin_dm, 0),
                    cap_strike: Dual::<1>::from_value(spec.cap_strike),
                    floor_strike: Dual::<1>::from_value(spec.floor_strike),
                    rate_normal_vol_sigma: Dual::<1>::from_value(spec.rate_normal_vol_sigma),
                    time_to_next_reset_tau1: Dual::<1>::from_value(spec.time_to_next_reset_tau1),
                },
            );
            let mut sp = spec;
            sp.discount_margin_dm += h;
            let mut sm = spec;
            sm.discount_margin_dm -= h;
            let fd = (evaluate_fabozzi_collared_frn_discount_margin(&sp).collared_frn_price
                - evaluate_fabozzi_collared_frn_discount_margin(&sm).collared_frn_price)
                / (2.0 * h);
            assert!((r_ad.collared_frn_price.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_fabozzi_collared_frn_discount_margin_module::*;

mod evaluate_tuckman_special_repo_fail_penalty_ctd_module {
    // Tuckman & Serrat (2011) & Fabozzi (2021) Special Repo Rate (`GC` vs. `Special`),
    // TMPG Fail-Penalty Floor & Cheapest-to-Deliver (`CTD`) Implied-Repo Switch Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct TuckmanSpecialRepoFailPenaltyCtdSpec<T: Real = f64> {
        pub general_collateral_rate_gc: T,
        pub tmpg_fail_penalty_rate_p: T,
        pub issue_short_demand_ratio_d: T,
        pub lendable_supply_elasticity_k: T,
        pub bond_clean_price: T,
        pub accrued_interest_at_spot: T,
        pub coupon_rate_annual: T,
        pub futures_invoice_price_at_delivery: T,
        pub holding_period_year_frac_tau: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct TuckmanSpecialRepoFailPenaltyCtdReport<T: Real = f64> {
        pub equilibrium_special_repo_rate: T,
        pub repo_specialness_spread_bps: T,
        pub tmpg_fail_floor_rate: T,
        pub financing_advantage_pv_per_100: T,
        pub implied_repo_rate_irr: T,
        pub net_basis_after_specialness: T,
        pub fail_probability_index: T,
    }

    pub fn evaluate_tuckman_special_repo_fail_penalty_ctd<T: Real>(
        spec: &TuckmanSpecialRepoFailPenaltyCtdSpec<T>,
    ) -> TuckmanSpecialRepoFailPenaltyCtdReport<T> {
        let gc = spec.general_collateral_rate_gc;
        let pen = spec.tmpg_fail_penalty_rate_p.max_val(T::from_f64(0.01));
        let tau = spec.holding_period_year_frac_tau.max_val(T::from_f64(0.01));

        // TMPG fail floor rate = -penalty (or GC - penalty in scarcity equilibrium):
        let fail_floor = gc - pen;
        let scarcity = (spec.issue_short_demand_ratio_d
            * spec.lendable_supply_elasticity_k.max_val(T::from_f64(0.01)))
        .max_val(T::zero());
        let fail_prob = T::one() - (-scarcity).exp();

        // Equilibrium special repo rate interpolates smoothly between GC and the TMPG fail floor:
        let special_rate = gc - pen * fail_prob;
        let specialness_bps = (gc - special_rate) * T::from_f64(1e4);

        let dirty_spot = (spec.bond_clean_price + spec.accrued_interest_at_spot).max_val(T::one());
        let fin_adv_pv = dirty_spot * (gc - special_rate) * tau;

        // Implied repo rate (IRR) from cash-and-carry to futures delivery:
        let carry_coupon = T::from_f64(100.0) * spec.coupon_rate_annual * tau;
        let proceeds = spec.futures_invoice_price_at_delivery + carry_coupon;
        let irr = ((proceeds - dirty_spot) / dirty_spot) / tau;

        // Net basis = DirtySpot * (1 + r_special * tau) - CarryCoupon - InvoicePrice:
        let net_basis = dirty_spot * (T::one() + special_rate * tau)
            - carry_coupon
            - spec.futures_invoice_price_at_delivery;

        TuckmanSpecialRepoFailPenaltyCtdReport {
            equilibrium_special_repo_rate: special_rate,
            repo_specialness_spread_bps: specialness_bps,
            tmpg_fail_floor_rate: fail_floor,
            financing_advantage_pv_per_100: fin_adv_pv,
            implied_repo_rate_irr: irr,
            net_basis_after_specialness: net_basis,
            fail_probability_index: fail_prob,
        }
    }

    #[cfg(test)]
    mod domain_suite_53_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn tuckman_special_repo_fail_penalty_ctd_invariants_and_ad_parity() {
            let spec = TuckmanSpecialRepoFailPenaltyCtdSpec {
                general_collateral_rate_gc: 0.048,
                tmpg_fail_penalty_rate_p: 0.030,
                issue_short_demand_ratio_d: 1.4,
                lendable_supply_elasticity_k: 0.55,
                bond_clean_price: 98.5,
                accrued_interest_at_spot: 0.85,
                coupon_rate_annual: 0.040,
                futures_invoice_price_at_delivery: 99.65,
                holding_period_year_frac_tau: 0.25,
            };
            let rep = evaluate_tuckman_special_repo_fail_penalty_ctd(&spec);
            assert!(rep.equilibrium_special_repo_rate < spec.general_collateral_rate_gc);
            assert!(rep.equilibrium_special_repo_rate > rep.tmpg_fail_floor_rate);
            assert!(rep.repo_specialness_spread_bps > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_tuckman_special_repo_fail_penalty_ctd(
                &TuckmanSpecialRepoFailPenaltyCtdSpec {
                    general_collateral_rate_gc: Dual::<1>::variable(
                        spec.general_collateral_rate_gc,
                        0,
                    ),
                    tmpg_fail_penalty_rate_p: Dual::<1>::from_value(spec.tmpg_fail_penalty_rate_p),
                    issue_short_demand_ratio_d: Dual::<1>::from_value(
                        spec.issue_short_demand_ratio_d,
                    ),
                    lendable_supply_elasticity_k: Dual::<1>::from_value(
                        spec.lendable_supply_elasticity_k,
                    ),
                    bond_clean_price: Dual::<1>::from_value(spec.bond_clean_price),
                    accrued_interest_at_spot: Dual::<1>::from_value(spec.accrued_interest_at_spot),
                    coupon_rate_annual: Dual::<1>::from_value(spec.coupon_rate_annual),
                    futures_invoice_price_at_delivery: Dual::<1>::from_value(
                        spec.futures_invoice_price_at_delivery,
                    ),
                    holding_period_year_frac_tau: Dual::<1>::from_value(
                        spec.holding_period_year_frac_tau,
                    ),
                },
            );
            let mut sp = spec;
            sp.general_collateral_rate_gc += h;
            let mut sm = spec;
            sm.general_collateral_rate_gc -= h;
            let fd = (evaluate_tuckman_special_repo_fail_penalty_ctd(&sp)
                .net_basis_after_specialness
                - evaluate_tuckman_special_repo_fail_penalty_ctd(&sm).net_basis_after_specialness)
                / (2.0 * h);
            assert!((r_ad.net_basis_after_specialness.d[0] - fd).abs() < 1e-6);
        }
    }
}
pub use evaluate_tuckman_special_repo_fail_penalty_ctd_module::*;

mod evaluate_fabozzi_callable_muni_ytw_refunding_module {
    // Fabozzi (2012) & Kalotay–Dorigan Tax-Exempt Callable Municipal Bond Yield-to-Worst
    // (`YTW`), Refunding Efficiency & Implied Marginal Tax-Rate (`Muni/UST`) Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct FabozziCallableMuniYtwRefundingSpec<T: Real = f64> {
        pub tax_exempt_coupon_rate_c: T,
        pub muni_mmd_benchmark_yield_ym: T,
        pub treasury_yield_y_ust: T,
        pub statutory_marginal_tax_rate_tau: T,
        pub call_protection_years_tc: T,
        pub final_maturity_years_t: T,
        pub call_price_par: T,
        pub muni_yield_volatility_sigma: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct FabozziCallableMuniYtwRefundingReport<T: Real = f64> {
        pub option_free_muni_bond_price: T,
        pub embedded_refunding_call_value: T,
        pub callable_muni_market_price: T,
        pub yield_to_call_ytc: T,
        pub yield_to_maturity_ytm: T,
        pub smooth_yield_to_worst_ytw: T,
        pub muni_to_treasury_yield_ratio: T,
        pub implied_breakeven_marginal_tax_rate: T,
        pub kalotay_refunding_efficiency: T,
    }

    fn bond_pv_and_dpv_cont<T: Real>(c: T, y: T, mat: T, par: T) -> (T, T) {
        let r = y.max_val(T::from_f64(1e-4));
        let df = (-r * mat).exp();
        let ann = (T::one() - df) / r;
        let pv = par * c * ann + par * df;
        let d_ann = (mat * df - ann) / r;
        let dpv_dy = par * c * d_ann - par * mat * df;
        (pv, dpv_dy)
    }

    fn solve_yield_cont<T: Real>(c: T, target_pv: T, mat: T, par: T, y_init: T) -> T {
        let mut y = y_init.max_val(T::from_f64(1e-3));
        for _ in 0..3 {
            let (pv, dpv) = bond_pv_and_dpv_cont(c, y, mat, par);
            let denom = dpv.min_val(T::from_f64(-1e-4));
            y = (y - (pv - target_pv) / denom).max_val(T::from_f64(1e-3));
        }
        y
    }

    pub fn evaluate_fabozzi_callable_muni_ytw_refunding<T: Real>(
        spec: &FabozziCallableMuniYtwRefundingSpec<T>,
    ) -> FabozziCallableMuniYtwRefundingReport<T> {
        let c = spec.tax_exempt_coupon_rate_c.max_val(T::from_f64(0.005));
        let ym = spec.muni_mmd_benchmark_yield_ym.max_val(T::from_f64(0.005));
        let y_ust = spec.treasury_yield_y_ust.max_val(T::from_f64(0.005));
        let tau_tax = spec
            .statutory_marginal_tax_rate_tau
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.60));
        let tc = spec.call_protection_years_tc.max_val(T::from_f64(0.5));
        let t_mat = spec.final_maturity_years_t.max_val(tc + T::from_f64(1.0));
        let par = spec.call_price_par.max_val(T::from_f64(50.0));
        let sig = spec.muni_yield_volatility_sigma.max_val(T::from_f64(0.02));

        // After-tax Treasury rate vs Muni rate wedge:
        let after_tax_ust = y_ust * (T::one() - tau_tax);
        let tax_exempt_disc = T::from_f64(0.75) * ym + T::from_f64(0.25) * after_tax_ust;

        let (p_opt_free, _) = bond_pv_and_dpv_cont(c, tax_exempt_disc, t_mat, par);

        // Forward bond value and exact forward modified duration at call date tc over remaining life (t_mat - tc):
        let rem_tenor = t_mat - tc;
        let (fwd_bond_tc, d_fwd_dy) = bond_pv_and_dpv_cont(c, tax_exempt_disc, rem_tenor, par);
        let fwd_mod_dur = (-d_fwd_dy) / fwd_bond_tc.max_val(T::from_f64(1e-4));
        let df_tc = (-tax_exempt_disc * tc).exp();

        // Black-76 swaption/bond-option on forward municipal bond at tc using exact duration mapping:
        let price_vol = sig * tax_exempt_disc * fwd_mod_dur;
        let std_b = (price_vol * tc.sqrt()).max_val(T::from_f64(1e-4));
        let d1 = ((fwd_bond_tc / par).ln() + T::from_f64(0.5) * std_b * std_b) / std_b;
        let d2 = d1 - std_b;
        let call_val = df_tc * (fwd_bond_tc * norm_cdf(d1) - par * norm_cdf(d2));
        let p_callable = p_opt_free - call_val;

        // Exact continuous-compounding YTC (over tc) and YTM (over t_mat) via 3 smooth Newton steps:
        let ann_tc = (T::one() - (-tax_exempt_disc * tc).exp()) / tax_exempt_disc;
        let ann_tm = (T::one() - (-tax_exempt_disc * t_mat).exp()) / tax_exempt_disc;
        let ytc_init = c - (p_callable - par) / (par * ann_tc);
        let ytm_init = c - (p_callable - par) / (par * ann_tm);
        let ytc = solve_yield_cont(c, p_callable, tc, par, ytc_init);
        let ytm = solve_yield_cont(c, p_callable, t_mat, par, ytm_init);

        let eps = T::from_f64(1e-4);
        let diff_y = ytc - ytm;
        let ytw = T::from_f64(0.5) * (ytc + ytm - (diff_y * diff_y + eps * eps).sqrt());

        let muni_ust_ratio = ytw / y_ust;
        let implied_tax_rate = T::one() - muni_ust_ratio;

        let intrinsic_ref_savings = df_tc
            * T::from_f64(0.5)
            * ((fwd_bond_tc - par)
                + ((fwd_bond_tc - par) * (fwd_bond_tc - par) + eps * eps).sqrt());
        let kalotay_eff = intrinsic_ref_savings / call_val.max_val(T::from_f64(1e-4));

        FabozziCallableMuniYtwRefundingReport {
            option_free_muni_bond_price: p_opt_free,
            embedded_refunding_call_value: call_val,
            callable_muni_market_price: p_callable,
            yield_to_call_ytc: ytc,
            yield_to_maturity_ytm: ytm,
            smooth_yield_to_worst_ytw: ytw,
            muni_to_treasury_yield_ratio: muni_ust_ratio,
            implied_breakeven_marginal_tax_rate: implied_tax_rate,
            kalotay_refunding_efficiency: kalotay_eff,
        }
    }

    #[cfg(test)]
    mod evaluate_fabozzi_callable_muni_ytw_refunding_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fabozzi_callable_muni_ytw_refunding_invariants_and_ad_parity() {
            let spec = FabozziCallableMuniYtwRefundingSpec {
                tax_exempt_coupon_rate_c: 0.050,
                muni_mmd_benchmark_yield_ym: 0.031,
                treasury_yield_y_ust: 0.044,
                statutory_marginal_tax_rate_tau: 0.37,
                call_protection_years_tc: 5.0,
                final_maturity_years_t: 10.0,
                call_price_par: 100.0,
                muni_yield_volatility_sigma: 0.18,
            };
            let rep = evaluate_fabozzi_callable_muni_ytw_refunding(&spec);
            assert!(rep.callable_muni_market_price < rep.option_free_muni_bond_price);
            assert!(rep.smooth_yield_to_worst_ytw <= rep.yield_to_maturity_ytm + 1e-4);
            assert!(
                rep.kalotay_refunding_efficiency > 0.5 && rep.kalotay_refunding_efficiency <= 1.0
            );

            let h = 1e-6;
            let r_ad = evaluate_fabozzi_callable_muni_ytw_refunding(
                &FabozziCallableMuniYtwRefundingSpec {
                    tax_exempt_coupon_rate_c: Dual::<8>::variable(spec.tax_exempt_coupon_rate_c, 0),
                    muni_mmd_benchmark_yield_ym: Dual::<8>::variable(
                        spec.muni_mmd_benchmark_yield_ym,
                        1,
                    ),
                    treasury_yield_y_ust: Dual::<8>::variable(spec.treasury_yield_y_ust, 2),
                    statutory_marginal_tax_rate_tau: Dual::<8>::variable(
                        spec.statutory_marginal_tax_rate_tau,
                        3,
                    ),
                    call_protection_years_tc: Dual::<8>::variable(spec.call_protection_years_tc, 4),
                    final_maturity_years_t: Dual::<8>::variable(spec.final_maturity_years_t, 5),
                    call_price_par: Dual::<8>::variable(spec.call_price_par, 6),
                    muni_yield_volatility_sigma: Dual::<8>::variable(
                        spec.muni_yield_volatility_sigma,
                        7,
                    ),
                },
            );
            for (idx, &ad_deriv) in r_ad.callable_muni_market_price.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.tax_exempt_coupon_rate_c += h;
                        sm.tax_exempt_coupon_rate_c -= h;
                    }
                    1 => {
                        sp.muni_mmd_benchmark_yield_ym += h;
                        sm.muni_mmd_benchmark_yield_ym -= h;
                    }
                    2 => {
                        sp.treasury_yield_y_ust += h;
                        sm.treasury_yield_y_ust -= h;
                    }
                    3 => {
                        sp.statutory_marginal_tax_rate_tau += h;
                        sm.statutory_marginal_tax_rate_tau -= h;
                    }
                    4 => {
                        sp.call_protection_years_tc += h;
                        sm.call_protection_years_tc -= h;
                    }
                    5 => {
                        sp.final_maturity_years_t += h;
                        sm.final_maturity_years_t -= h;
                    }
                    6 => {
                        sp.call_price_par += h;
                        sm.call_price_par -= h;
                    }
                    _ => {
                        sp.muni_yield_volatility_sigma += h;
                        sm.muni_yield_volatility_sigma -= h;
                    }
                }
                let fd = (evaluate_fabozzi_callable_muni_ytw_refunding(&sp)
                    .callable_muni_market_price
                    - evaluate_fabozzi_callable_muni_ytw_refunding(&sm).callable_muni_market_price)
                    / (2.0 * h);
                assert!((ad_deriv - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_fabozzi_callable_muni_ytw_refunding_module::*;

mod evaluate_tuckman_par_spot_forward_rolldown_horizon_return_module {
    // Tuckman & Serrat (2012, Fixed Income Securities, 3rd ed., Ch. 4-5) & Fabozzi (2021)
    // Par-Spot-Forward Roll-Down Horizon Return, Carry-Convexity Pickup & Breakeven Yield Widening Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct TuckmanRolldownHorizonReturnSpec<T: Real = f64> {
        pub initial_maturity_years_t: T,
        pub holding_horizon_years_dt: T,
        pub current_bond_yield_y_t: T,
        pub rolled_maturity_yield_y_t_minus_dt: T,
        pub short_financing_rate_r_dt: T,
        pub annual_coupon_rate_c: T,
        pub modified_duration_at_horizon_d_mod: T,
        pub convexity_at_horizon_conv: T,
        pub expected_parallel_shift_bps: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct TuckmanRolldownHorizonReturnReport<T: Real = f64> {
        pub implied_forward_rate_f_dt_t: T,
        pub annualized_rolldown_yield_pickup_bps: T,
        pub rolldown_price_appreciation_pct: T,
        pub convexity_return_pickup_pct: T,
        pub horizon_expected_excess_return_pct: T,
        pub breakeven_yield_widening_bps: T,
    }

    pub fn evaluate_tuckman_par_spot_forward_rolldown_horizon_return<T: Real>(
        spec: &TuckmanRolldownHorizonReturnSpec<T>,
    ) -> TuckmanRolldownHorizonReturnReport<T> {
        let t = spec.initial_maturity_years_t.max_val(T::from_f64(1.0));
        let dt = spec
            .holding_horizon_years_dt
            .max_val(T::from_f64(0.05))
            .min_val(t - T::from_f64(0.25));
        let y_t = spec.current_bond_yield_y_t;
        let y_roll = spec.rolled_maturity_yield_y_t_minus_dt;
        let r_fin = spec.short_financing_rate_r_dt;
        let c = spec.annual_coupon_rate_c;
        let d_mod = spec
            .modified_duration_at_horizon_d_mod
            .max_val(T::from_f64(0.25));
        let conv = spec.convexity_at_horizon_conv.max_val(T::zero());
        let dy_view = spec.expected_parallel_shift_bps * T::from_f64(1e-4);

        let rem_t = (t - dt).max_val(T::from_f64(0.25));
        let f_impl = (y_t * t - r_fin * dt) / rem_t;

        let rolldown_dy = y_t - y_roll;
        let rolldown_bps = rolldown_dy * T::from_f64(10_000.0);
        let rolldown_price_pct = d_mod * rolldown_dy * T::from_f64(100.0);

        let dy_total = (y_roll - y_t) + dy_view;
        let conv_pct = T::from_f64(0.5) * conv * dy_total * dy_total * T::from_f64(100.0);

        let y_safe = y_t.max_val(T::from_f64(1e-4));
        let ann_t = (T::one() - (-y_safe * t).exp()) / y_safe;
        let ann_rem = (T::one() - (-y_safe * rem_t).exp()) / y_safe;
        let p_0 = T::one() + (c - y_t) * ann_t;
        let p_const_y = T::one() + (c - y_t) * ann_rem;
        let carry_excess = (c * dt + (p_const_y - p_0)) / p_0 - r_fin * dt;
        let total_excess_pct = (carry_excess - d_mod * dy_total
            + T::from_f64(0.5) * conv * dy_total * dy_total)
            * T::from_f64(100.0);

        let breakeven_bps = ((carry_excess + d_mod * rolldown_dy) / d_mod) * T::from_f64(10_000.0);

        TuckmanRolldownHorizonReturnReport {
            implied_forward_rate_f_dt_t: f_impl,
            annualized_rolldown_yield_pickup_bps: rolldown_bps,
            rolldown_price_appreciation_pct: rolldown_price_pct,
            convexity_return_pickup_pct: conv_pct,
            horizon_expected_excess_return_pct: total_excess_pct,
            breakeven_yield_widening_bps: breakeven_bps,
        }
    }

    #[cfg(test)]
    mod evaluate_tuckman_par_spot_forward_rolldown_horizon_return_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn tuckman_par_spot_forward_rolldown_invariants_and_ad_parity() {
            let spec = TuckmanRolldownHorizonReturnSpec {
                initial_maturity_years_t: 7.0,
                holding_horizon_years_dt: 1.0,
                current_bond_yield_y_t: 0.0445,
                rolled_maturity_yield_y_t_minus_dt: 0.0420,
                short_financing_rate_r_dt: 0.0380,
                annual_coupon_rate_c: 0.0450,
                modified_duration_at_horizon_d_mod: 5.20,
                convexity_at_horizon_conv: 32.5,
                expected_parallel_shift_bps: 15.0,
            };
            let rep = evaluate_tuckman_par_spot_forward_rolldown_horizon_return(&spec);
            assert!(rep.annualized_rolldown_yield_pickup_bps > 0.0);
            assert!(rep.rolldown_price_appreciation_pct > 0.0);
            assert!(rep.breakeven_yield_widening_bps > 20.0);

            let h = 1e-5;
            let r_ad = evaluate_tuckman_par_spot_forward_rolldown_horizon_return(
                &TuckmanRolldownHorizonReturnSpec {
                    initial_maturity_years_t: Dual::<9>::variable(spec.initial_maturity_years_t, 0),
                    holding_horizon_years_dt: Dual::<9>::variable(spec.holding_horizon_years_dt, 1),
                    current_bond_yield_y_t: Dual::<9>::variable(spec.current_bond_yield_y_t, 2),
                    rolled_maturity_yield_y_t_minus_dt: Dual::<9>::variable(
                        spec.rolled_maturity_yield_y_t_minus_dt,
                        3,
                    ),
                    short_financing_rate_r_dt: Dual::<9>::variable(
                        spec.short_financing_rate_r_dt,
                        4,
                    ),
                    annual_coupon_rate_c: Dual::<9>::variable(spec.annual_coupon_rate_c, 5),
                    modified_duration_at_horizon_d_mod: Dual::<9>::variable(
                        spec.modified_duration_at_horizon_d_mod,
                        6,
                    ),
                    convexity_at_horizon_conv: Dual::<9>::variable(
                        spec.convexity_at_horizon_conv,
                        7,
                    ),
                    expected_parallel_shift_bps: Dual::<9>::variable(
                        spec.expected_parallel_shift_bps,
                        8,
                    ),
                },
            );
            for (idx, &ad) in r_ad.horizon_expected_excess_return_pct.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.initial_maturity_years_t += h;
                        sm.initial_maturity_years_t -= h;
                    }
                    1 => {
                        sp.holding_horizon_years_dt += h;
                        sm.holding_horizon_years_dt -= h;
                    }
                    2 => {
                        sp.current_bond_yield_y_t += h;
                        sm.current_bond_yield_y_t -= h;
                    }
                    3 => {
                        sp.rolled_maturity_yield_y_t_minus_dt += h;
                        sm.rolled_maturity_yield_y_t_minus_dt -= h;
                    }
                    4 => {
                        sp.short_financing_rate_r_dt += h;
                        sm.short_financing_rate_r_dt -= h;
                    }
                    5 => {
                        sp.annual_coupon_rate_c += h;
                        sm.annual_coupon_rate_c -= h;
                    }
                    6 => {
                        sp.modified_duration_at_horizon_d_mod += h;
                        sm.modified_duration_at_horizon_d_mod -= h;
                    }
                    7 => {
                        sp.convexity_at_horizon_conv += h;
                        sm.convexity_at_horizon_conv -= h;
                    }
                    _ => {
                        sp.expected_parallel_shift_bps += h;
                        sm.expected_parallel_shift_bps -= h;
                    }
                }
                let rp = evaluate_tuckman_par_spot_forward_rolldown_horizon_return(&sp);
                let rm = evaluate_tuckman_par_spot_forward_rolldown_horizon_return(&sm);
                let fd = (rp.horizon_expected_excess_return_pct
                    - rm.horizon_expected_excess_return_pct)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_tuckman_par_spot_forward_rolldown_horizon_return_module::*;

mod evaluate_kalotay_fabozzi_callable_bond_refunding_efficiency_module {
    // Kalotay, Williams & Fabozzi (1993, Financial Analysts Journal 49(3), "A Model for
    // Valuing Bonds and Embedded Options"), Kalotay, Yang & Fabozzi (2007, J. Fixed Income),
    // & Fabozzi (2021, The Handbook of Fixed Income Securities, 9th ed., Ch. 16 & 40)
    // Callable Bond Refunding NPV, Forfeited Call Time-Value & Refunding Efficiency Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct KalotayFabozziRefundingEfficiencySpec<T: Real = f64> {
        pub notional_par_n: T,
        pub existing_coupon_rate_c_old: T,
        pub refunding_market_yield_y_ref: T,
        pub call_strike_price_k_call: T,
        pub refunding_issuance_cost_c_iss: T,
        pub yield_volatility_sigma_y: T,
        pub deferral_option_horizon_t_c: T,
        pub remaining_maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct KalotayFabozziRefundingEfficiencyReport<T: Real = f64> {
        pub pv_gross_interest_savings: T,
        pub net_present_value_refunding_savings: T,
        pub forfeited_waiting_option_time_value: T,
        pub total_embedded_refunding_call_value: T,
        pub kalotay_refunding_efficiency_ratio: T,
        pub breakeven_refunding_yield_rate: T,
    }

    pub fn evaluate_kalotay_fabozzi_callable_bond_refunding_efficiency<T: Real>(
        spec: &KalotayFabozziRefundingEfficiencySpec<T>,
    ) -> KalotayFabozziRefundingEfficiencyReport<T> {
        let par = spec.notional_par_n.max_val(T::from_f64(1.0));
        let c_old = spec.existing_coupon_rate_c_old.max_val(T::from_f64(0.005));
        let y_ref = spec
            .refunding_market_yield_y_ref
            .max_val(T::from_f64(0.002));
        let k_call = spec.call_strike_price_k_call.max_val(par);
        let c_iss = spec.refunding_issuance_cost_c_iss.max_val(T::from_f64(0.0));
        let sig_y = spec.yield_volatility_sigma_y.max_val(T::from_f64(0.02));
        let t_c = spec.deferral_option_horizon_t_c.max_val(T::from_f64(0.25));
        let t_mat = spec
            .remaining_maturity_years_t
            .max_val(t_c + T::from_f64(1.0));

        let annuity = (T::one() - (-y_ref * t_mat).exp()) / y_ref;
        let gross_savings = par * (c_old - y_ref) * annuity;
        let call_prem_and_cost = (k_call - par) + c_iss;
        let npv_savings = gross_savings - call_prem_and_cost;

        // Black-76 receiver swaption / bond call value over deferral window t_c on remaining tail (t_mat - t_c):
        let tail_t = t_mat - t_c;
        let tail_ann = (T::one() - (-y_ref * tail_t).exp()) / y_ref;
        let k_yield = (c_old - call_prem_and_cost / (par * tail_ann)).max_val(T::from_f64(0.002));
        let sqrt_tc = t_c.sqrt();
        let std_y = sig_y * sqrt_tc;
        let d1 = ((k_yield / y_ref).ln() + T::from_f64(0.5) * sig_y * sig_y * t_c) / std_y;
        let d2 = d1 - std_y;

        let df_tc = (-y_ref * t_c).exp();
        // By Black-76 put-call parity, the pure extrinsic waiting time value over t_c on the
        // remaining tail (t_mat - t_c) is strictly positive:
        // V_rec(k_yield, y_ref) - Intrinsic_fwd = df_tc * Par * A_tail * (y_ref * N(-d2) - k_yield * N(-d1)) > 0
        let time_val = df_tc * par * tail_ann * (y_ref * norm_cdf(-d2) - k_yield * norm_cdf(-d1));
        let total_opt_val = npv_savings + time_val;

        let eff_ratio = npv_savings / total_opt_val.max_val(T::from_f64(1e-4));
        let breakeven_y = c_old - call_prem_and_cost / (par * annuity);

        KalotayFabozziRefundingEfficiencyReport {
            pv_gross_interest_savings: gross_savings,
            net_present_value_refunding_savings: npv_savings,
            forfeited_waiting_option_time_value: time_val,
            total_embedded_refunding_call_value: total_opt_val,
            kalotay_refunding_efficiency_ratio: eff_ratio,
            breakeven_refunding_yield_rate: breakeven_y,
        }
    }

    #[cfg(test)]
    mod evaluate_kalotay_fabozzi_callable_bond_refunding_efficiency_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn kalotay_fabozzi_refunding_efficiency_invariants_and_ad_parity() {
            let spec = KalotayFabozziRefundingEfficiencySpec {
                notional_par_n: 100.0,
                existing_coupon_rate_c_old: 0.055,
                refunding_market_yield_y_ref: 0.038,
                call_strike_price_k_call: 101.0,
                refunding_issuance_cost_c_iss: 0.65,
                yield_volatility_sigma_y: 0.22,
                deferral_option_horizon_t_c: 1.50,
                remaining_maturity_years_t: 10.0,
            };
            let rep = evaluate_kalotay_fabozzi_callable_bond_refunding_efficiency(&spec);
            assert!(rep.net_present_value_refunding_savings > 0.0);
            assert!(
                rep.kalotay_refunding_efficiency_ratio > 0.50
                    && rep.kalotay_refunding_efficiency_ratio < 1.0
            );

            let h = 1e-6;
            let r_ad = evaluate_kalotay_fabozzi_callable_bond_refunding_efficiency(
                &KalotayFabozziRefundingEfficiencySpec {
                    notional_par_n: Dual::from_f64(spec.notional_par_n),
                    existing_coupon_rate_c_old: Dual::<6>::variable(
                        spec.existing_coupon_rate_c_old,
                        0,
                    ),
                    refunding_market_yield_y_ref: Dual::<6>::variable(
                        spec.refunding_market_yield_y_ref,
                        1,
                    ),
                    call_strike_price_k_call: Dual::<6>::variable(spec.call_strike_price_k_call, 2),
                    refunding_issuance_cost_c_iss: Dual::<6>::variable(
                        spec.refunding_issuance_cost_c_iss,
                        3,
                    ),
                    yield_volatility_sigma_y: Dual::<6>::variable(spec.yield_volatility_sigma_y, 4),
                    deferral_option_horizon_t_c: Dual::<6>::variable(
                        spec.deferral_option_horizon_t_c,
                        5,
                    ),
                    remaining_maturity_years_t: Dual::from_f64(spec.remaining_maturity_years_t),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.existing_coupon_rate_c_old,
                    spec.refunding_market_yield_y_ref,
                    spec.call_strike_price_k_call,
                    spec.refunding_issuance_cost_c_iss,
                    spec.yield_volatility_sigma_y,
                    spec.deferral_option_horizon_t_c,
                ];
                v[i] += h;
                let up = evaluate_kalotay_fabozzi_callable_bond_refunding_efficiency(
                    &KalotayFabozziRefundingEfficiencySpec {
                        existing_coupon_rate_c_old: v[0],
                        refunding_market_yield_y_ref: v[1],
                        call_strike_price_k_call: v[2],
                        refunding_issuance_cost_c_iss: v[3],
                        yield_volatility_sigma_y: v[4],
                        deferral_option_horizon_t_c: v[5],
                        ..spec
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_kalotay_fabozzi_callable_bond_refunding_efficiency(
                    &KalotayFabozziRefundingEfficiencySpec {
                        existing_coupon_rate_c_old: v[0],
                        refunding_market_yield_y_ref: v[1],
                        call_strike_price_k_call: v[2],
                        refunding_issuance_cost_c_iss: v[3],
                        yield_volatility_sigma_y: v[4],
                        deferral_option_horizon_t_c: v[5],
                        ..spec
                    },
                );
                let fd = (up.kalotay_refunding_efficiency_ratio
                    - dn.kalotay_refunding_efficiency_ratio)
                    / (2.0 * h);
                assert!(
                    (r_ad.kalotay_refunding_efficiency_ratio.d[i] - fd).abs() < 1e-5,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_kalotay_fabozzi_callable_bond_refunding_efficiency_module::*;

pub mod black_derman_toy_callable_suite65_mod {
    use ql_ad::Real;

    /// Specification for the Black–Derman–Toy (1990) lognormal short-rate lattice callable bond pricer.
    ///
    /// Reference: Black, Derman & Toy (1990, *Financial Analysts Journal* 46(1), "A One-Factor Model
    /// of Interest Rates and Its Application to Treasury Bond Options"), Tuckman & Serrat (2022, Ch. 11).
    #[derive(Debug, Clone, Copy)]
    pub struct BlackDermanToyCallableBondSpec<T: Real> {
        /// Bond par value $F_{\text{par}}$.
        pub par_value: T,
        /// Annualized coupon rate $c$.
        pub coupon_rate_c: T,
        /// Initial short-rate median level $r_0$.
        pub initial_short_rate_r0: T,
        /// Term structure drift slope $s_{\text{slope}}$ per year.
        pub term_slope_s: T,
        /// Lognormal short-rate volatility $\sigma_r$.
        pub short_rate_log_vol_sigma: T,
        /// Option-Adjusted Spread (OAS) $s_{\text{oas}}$.
        pub option_adjusted_spread_oas: T,
        /// Issuer call price $K_{\text{call}}$.
        pub call_price_k: T,
        /// Maturity $T$ (in years).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_black_derman_toy_lognormal_short_rate_callable_bond`].
    #[derive(Debug, Clone, Copy)]
    pub struct BlackDermanToyCallableBondAnalytics<T: Real> {
        /// Option-free straight bond price on the BDT lattice.
        pub straight_bond_price: T,
        /// Callable bond price on the BDT lattice with smooth call exercise.
        pub bdt_callable_bond_price: T,
        /// Embedded issuer call option value $V_{\text{call}} = B_{\text{straight}} - B_{\text{callable}}$.
        pub embedded_call_option_value: T,
        /// Effective BDT lattice rate sensitivity (modified duration proxy).
        pub effective_callable_duration: T,
    }

    /// Evaluates an option-free and callable bond on a 4-step recombining Black–Derman–Toy (1990)
    /// lognormal short-rate binomial tree with $N = 4$ discounted period coupons.
    pub fn evaluate_black_derman_toy_lognormal_short_rate_callable_bond<T: Real>(
        spec: &BlackDermanToyCallableBondSpec<T>,
    ) -> BlackDermanToyCallableBondAnalytics<T> {
        const STEPS: usize = 4;
        let t_mat = spec.maturity_years_t.max_val(T::from_f64(1e-4));
        let dt = t_mat / T::from_f64(STEPS as f64);
        let sqrt_dt = dt.sqrt();
        let sig = spec.short_rate_log_vol_sigma.max_val(T::from_f64(1e-4));
        let step_coupon = spec.par_value * spec.coupon_rate_c * dt;

        let smooth_pos = |x: T| -> T {
            let eps = spec.par_value * T::from_f64(0.025);
            T::from_f64(0.5) * (x + (x * x + eps * eps).sqrt())
        };

        let eval_tree = |rate_shift: T| -> (T, T) {
            let mut straight = [spec.par_value; STEPS + 1];
            let mut callable = [spec.par_value; STEPS + 1];

            for n in (0..STEPS).rev() {
                let t_n = T::from_f64(n as f64) * dt;
                let median_r = (spec.initial_short_rate_r0 + spec.term_slope_s * t_n + rate_shift)
                    .max_val(T::from_f64(1e-4));
                for j in 0..=n {
                    let node_shift = (2 * j as isize - n as isize) as f64;
                    let r_nj = median_r * (T::from_f64(node_shift) * sig * sqrt_dt).exp();
                    let df_nj = (-(r_nj + spec.option_adjusted_spread_oas) * dt).exp();

                    let s_cont =
                        df_nj * (T::from_f64(0.5) * (straight[j + 1] + straight[j]) + step_coupon);
                    straight[j] = s_cont;

                    let c_cont =
                        df_nj * (T::from_f64(0.5) * (callable[j + 1] + callable[j]) + step_coupon);
                    callable[j] = c_cont - smooth_pos(c_cont - spec.call_price_k);
                }
            }
            (straight[0], callable[0])
        };

        let (s0, c0) = eval_tree(T::zero());
        let dr = T::from_f64(0.001);
        let (_s_up, c_up) = eval_tree(dr);
        let (_s_dn, c_dn) = eval_tree(-dr);
        let eff_dur = (c_dn - c_up) / (T::from_f64(2.0) * dr * c0.max_val(T::from_f64(1.0)));

        BlackDermanToyCallableBondAnalytics {
            straight_bond_price: s0,
            bdt_callable_bond_price: c0,
            embedded_call_option_value: s0 - c0,
            effective_callable_duration: eff_dur,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bdt_callable_bond_invariants_and_ad_parity() {
            let base = BlackDermanToyCallableBondSpec {
                par_value: 100.0_f64,
                coupon_rate_c: 0.055,
                initial_short_rate_r0: 0.038,
                term_slope_s: 0.002,
                short_rate_log_vol_sigma: 0.22,
                option_adjusted_spread_oas: 0.006,
                call_price_k: 101.0,
                maturity_years_t: 4.0,
            };
            let r0 = evaluate_black_derman_toy_lognormal_short_rate_callable_bond(&base);
            assert!(r0.straight_bond_price > r0.bdt_callable_bond_price);
            assert!(r0.embedded_call_option_value > 0.0);
            assert!(r0.effective_callable_duration > 0.5);

            let inputs = [
                base.par_value,
                base.coupon_rate_c,
                base.initial_short_rate_r0,
                base.term_slope_s,
                base.short_rate_log_vol_sigma,
                base.option_adjusted_spread_oas,
                base.call_price_k,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BlackDermanToyCallableBondSpec {
                par_value: dual_in[0],
                coupon_rate_c: dual_in[1],
                initial_short_rate_r0: dual_in[2],
                term_slope_s: dual_in[3],
                short_rate_log_vol_sigma: dual_in[4],
                option_adjusted_spread_oas: dual_in[5],
                call_price_k: dual_in[6],
                maturity_years_t: dual_in[7],
            };
            let r_ad = evaluate_black_derman_toy_lognormal_short_rate_callable_bond(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| BlackDermanToyCallableBondSpec {
                    par_value: v[0],
                    coupon_rate_c: v[1],
                    initial_short_rate_r0: v[2],
                    term_slope_s: v[3],
                    short_rate_log_vol_sigma: v[4],
                    option_adjusted_spread_oas: v[5],
                    call_price_k: v[6],
                    maturity_years_t: v[7],
                };
                let rp = evaluate_black_derman_toy_lognormal_short_rate_callable_bond(&mk(p));
                let rm = evaluate_black_derman_toy_lognormal_short_rate_callable_bond(&mk(m));
                let fd = (rp.bdt_callable_bond_price - rm.bdt_callable_bond_price) / (2.0 * eps);
                assert!(
                    r_ad.bdt_callable_bond_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.bdt_callable_bond_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use black_derman_toy_callable_suite65_mod::*;

pub mod tuckman_corporate_make_whole_bond_mod {
    use ql_ad::Real;
    use ql_math::{norm_cdf, norm_pdf};

    /// Specification for the Tuckman & Serrat (2022, Ch. 19) / Fabozzi (2021, Ch. 16)
    /// corporate bond with a make-whole call provision and a late-stage par-call window.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct TuckmanMakeWholeCallableBondSpec<T: Real> {
        /// Bond par face value $F > 0$.
        pub par_value_f: T,
        /// Annual contractual coupon rate $c > 0$.
        pub coupon_rate_c: T,
        /// Benchmark risk-free Treasury yield $y_{\text{tsy}}$.
        pub treasury_benchmark_yield_y: T,
        /// Issuer corporate credit spread $s_{\text{corp}} > 0$ over Treasury.
        pub issuer_credit_spread_s_corp: T,
        /// Contractual make-whole discount spread $s_{\text{mw}} > 0$ over Treasury.
        pub make_whole_treasury_spread_s_mw: T,
        /// Normal volatility $\sigma_s > 0$ of the issuer's credit spread.
        pub credit_spread_volatility_sigma_s: T,
        /// Normal volatility $\sigma_y > 0$ of the benchmark Treasury yield.
        pub treasury_yield_volatility_sigma_y: T,
        /// Make-whole evaluation horizon $T_{\text{mw}} > 0$ (years).
        pub make_whole_horizon_t_mw: T,
        /// Par-call window start date $T_{\text{par}} > T_{\text{mw}}$ (years).
        pub par_call_window_start_t_par: T,
        /// Final stated bond maturity $T_{\text{mat}} > T_{\text{par}}$ (years).
        pub maturity_years_t_mat: T,
    }

    /// Analytical outputs for the Tuckman corporate make-whole & par-call bond engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct TuckmanMakeWholeCallableBondAnalytics<T: Real> {
        /// Straight non-callable corporate bond present value.
        pub straight_corporate_bond_price: T,
        /// Contractual make-whole redemption price at current Treasury yield + $s_{\text{mw}}$.
        pub make_whole_redemption_price: T,
        /// Combined value of the embedded make-whole spread option and late par-call option.
        pub embedded_make_whole_and_par_call_value: T,
        /// Callable corporate bond price (straight bond minus embedded issuer call options).
        pub callable_corporate_bond_price: T,
    }

    fn risky_bond_pv<T: Real>(par: T, c: T, yield_rate: T, tau: T) -> T {
        let y = yield_rate.max_val(T::from_f64(1e-4));
        let df = (-y * tau).exp();
        let ann = (T::one() - df) / y;
        par * c * ann + par * df
    }

    /// Evaluates a corporate bond with a make-whole call at $T_{\text{mw}}$ (driven by credit-spread
    /// compression below $s_{\text{mw}}$) and a par-call window at $T_{\text{par}}$.
    pub fn evaluate_tuckman_corporate_make_whole_and_par_call_bond<T: Real>(
        spec: &TuckmanMakeWholeCallableBondSpec<T>,
    ) -> TuckmanMakeWholeCallableBondAnalytics<T> {
        let par = spec.par_value_f.max_val(T::from_f64(1.0));
        let t_mat = spec.maturity_years_t_mat.max_val(T::from_f64(0.5));
        let t_par = spec
            .par_call_window_start_t_par
            .max_val(T::from_f64(0.25))
            .min_val(t_mat - T::from_f64(0.1));
        let t_mw = spec
            .make_whole_horizon_t_mw
            .max_val(T::from_f64(0.1))
            .min_val(t_par - T::from_f64(0.1));

        let y_corp = spec.treasury_benchmark_yield_y + spec.issuer_credit_spread_s_corp;
        let y_mw = spec.treasury_benchmark_yield_y + spec.make_whole_treasury_spread_s_mw;

        let straight_pv = risky_bond_pv(par, spec.coupon_rate_c, y_corp, t_mat);
        let rem_mw = t_mat - t_mw;
        let mw_redemption = risky_bond_pv(par, spec.coupon_rate_c, y_mw, rem_mw).max_val(par);

        // 1. Make-whole call at T_mw: exercised only if issuer credit spread s_corp(T_mw) compresses below s_mw.
        // Forward bond value at T_mw is sensitivity-scaled by risky duration D_mw * par.
        let dur_mw = (T::one() - (-y_corp * rem_mw).exp()) / y_corp.max_val(T::from_f64(1e-4));
        let sig_s = spec
            .credit_spread_volatility_sigma_s
            .max_val(T::from_f64(1e-5));
        let std_s_mw = sig_s * t_mw.sqrt();
        let spread_diff_mw =
            spec.make_whole_treasury_spread_s_mw - spec.issuer_credit_spread_s_corp;
        let d_mw = spread_diff_mw / std_s_mw;
        let df_mw = (-y_corp * t_mw).exp();
        let mw_call_val =
            df_mw * par * dur_mw * (spread_diff_mw * norm_cdf(d_mw) + std_s_mw * norm_pdf(d_mw));

        // 2. Par-call window at T_par (e.g. 3-6 months prior to maturity): standard call at par F
        // driven by total corporate yield volatility sqrt(sigma_y^2 + sigma_s^2).
        let rem_par = t_mat - t_par;
        let dur_par = (T::one() - (-y_corp * rem_par).exp()) / y_corp.max_val(T::from_f64(1e-4));
        let sig_y = spec
            .treasury_yield_volatility_sigma_y
            .max_val(T::from_f64(1e-5));
        let sig_tot = (sig_y * sig_y + sig_s * sig_s).sqrt();
        let std_y_par = sig_tot * t_par.sqrt();
        let yield_in_the_money = spec.coupon_rate_c - y_corp;
        let d_par = yield_in_the_money / std_y_par;
        let df_par = (-y_corp * t_par).exp();
        let par_call_val = df_par
            * par
            * dur_par
            * (yield_in_the_money * norm_cdf(d_par) + std_y_par * norm_pdf(d_par));

        let total_call = mw_call_val + par_call_val;

        TuckmanMakeWholeCallableBondAnalytics {
            straight_corporate_bond_price: straight_pv,
            make_whole_redemption_price: mw_redemption,
            embedded_make_whole_and_par_call_value: total_call,
            callable_corporate_bond_price: straight_pv - total_call,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn tuckman_make_whole_callable_bond_invariants_and_ad_parity() {
            let base = TuckmanMakeWholeCallableBondSpec {
                par_value_f: 100.0_f64,
                coupon_rate_c: 0.056,
                treasury_benchmark_yield_y: 0.038,
                issuer_credit_spread_s_corp: 0.014,
                make_whole_treasury_spread_s_mw: 0.0030,
                credit_spread_volatility_sigma_s: 0.0065,
                treasury_yield_volatility_sigma_y: 0.0085,
                make_whole_horizon_t_mw: 2.0,
                par_call_window_start_t_par: 4.5,
                maturity_years_t_mat: 5.0,
            };
            let r0 = evaluate_tuckman_corporate_make_whole_and_par_call_bond(&base);
            assert!(r0.make_whole_redemption_price > r0.straight_corporate_bond_price);
            assert!(r0.embedded_make_whole_and_par_call_value > 0.0);
            assert!(r0.callable_corporate_bond_price < r0.straight_corporate_bond_price);

            let inputs = [
                base.par_value_f,
                base.coupon_rate_c,
                base.treasury_benchmark_yield_y,
                base.issuer_credit_spread_s_corp,
                base.make_whole_treasury_spread_s_mw,
                base.credit_spread_volatility_sigma_s,
                base.treasury_yield_volatility_sigma_y,
                base.make_whole_horizon_t_mw,
                base.par_call_window_start_t_par,
                base.maturity_years_t_mat,
            ];
            let mut dual_in = [Dual::<10>::from_f64(0.0); 10];
            for i in 0..10 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = TuckmanMakeWholeCallableBondSpec {
                par_value_f: dual_in[0],
                coupon_rate_c: dual_in[1],
                treasury_benchmark_yield_y: dual_in[2],
                issuer_credit_spread_s_corp: dual_in[3],
                make_whole_treasury_spread_s_mw: dual_in[4],
                credit_spread_volatility_sigma_s: dual_in[5],
                treasury_yield_volatility_sigma_y: dual_in[6],
                make_whole_horizon_t_mw: dual_in[7],
                par_call_window_start_t_par: dual_in[8],
                maturity_years_t_mat: dual_in[9],
            };
            let r_ad = evaluate_tuckman_corporate_make_whole_and_par_call_bond(&spec_ad);
            let eps = 1e-6;
            for i in 0..10 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 10]| TuckmanMakeWholeCallableBondSpec {
                    par_value_f: v[0],
                    coupon_rate_c: v[1],
                    treasury_benchmark_yield_y: v[2],
                    issuer_credit_spread_s_corp: v[3],
                    make_whole_treasury_spread_s_mw: v[4],
                    credit_spread_volatility_sigma_s: v[5],
                    treasury_yield_volatility_sigma_y: v[6],
                    make_whole_horizon_t_mw: v[7],
                    par_call_window_start_t_par: v[8],
                    maturity_years_t_mat: v[9],
                };
                let rp = evaluate_tuckman_corporate_make_whole_and_par_call_bond(&mk(p));
                let rm = evaluate_tuckman_corporate_make_whole_and_par_call_bond(&mk(m));
                let fd = (rp.callable_corporate_bond_price - rm.callable_corporate_bond_price)
                    / (2.0 * eps);
                assert!(
                    r_ad.callable_corporate_bond_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.callable_corporate_bond_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use tuckman_corporate_make_whole_bond_mod::*;
mod kalotay_williams_fabozzi_callable_mod {
    use ql_ad::Real;

    /// Input specification for the Kalotay, Williams & Fabozzi (1993, *FAJ*) lognormal
    /// short-rate lattice callable bond valuation and refunding efficiency model.
    #[derive(Debug, Clone, Copy)]
    pub struct KalotayWilliamsFabozziCallableSpec<T: Real> {
        /// Initial 1-period benchmark short rate $r_0$ (e.g., $0.042$).
        pub initial_short_rate_r0: T,
        /// Annualized short-rate drift $\mu_r$ (e.g., $0.005$).
        pub short_rate_drift_mu: T,
        /// Lognormal short-rate volatility $\sigma_r$ (e.g., $0.18$).
        pub lognormal_short_rate_volatility_sigma: T,
        /// Option-Adjusted Spread (OAS) $s_{\text{OAS}}$ over benchmark short rates (e.g., $0.011$).
        pub option_adjusted_spread_oas: T,
        /// Annual coupon rate $c$ on par $100$ face value (e.g., $0.058$).
        pub coupon_rate_c: T,
        /// Call strike price $K_{\text{call}}$ per $100$ par (e.g., $100.5$).
        pub call_strike_price_k: T,
        /// Refunding issuance transaction cost per $100$ par $c_{\text{ref}}$ (e.g., $0.65$).
        pub refunding_transaction_cost: T,
        /// Lattice period step $\Delta t$ in years (e.g., $1.0$).
        pub step_years_dt: T,
    }

    /// Analytics returned by [`evaluate_kalotay_williams_fabozzi_callable_bond_refunding_efficiency`].
    #[derive(Debug, Clone, Copy)]
    pub struct KalotayWilliamsFabozziCallableAnalytics<T: Real> {
        /// Straight (non-callable) option-free bond value per $100$ par.
        pub straight_option_free_bond_value: T,
        /// Callable bond value per $100$ par under OAS discounting.
        pub callable_bond_value: T,
        /// Embedded issuer call option value $V_{\text{str}} - V_{\text{call}}$.
        pub embedded_call_option_value: T,
        /// Kalotay–Williams–Fabozzi net refunding efficiency ratio $\text{PV}(\text{Net Savings}) / V_{\text{opt}}$.
        pub refunding_efficiency_ratio: T,
    }

    /// Evaluates a 2-period Kalotay, Williams & Fabozzi (1993) lognormal short-rate lattice
    /// for a callable corporate/municipal bond, returning the option-free value, callable value,
    /// embedded call option value, and net refunding efficiency.
    pub fn evaluate_kalotay_williams_fabozzi_callable_bond_refunding_efficiency<T: Real>(
        spec: &KalotayWilliamsFabozziCallableSpec<T>,
    ) -> KalotayWilliamsFabozziCallableAnalytics<T> {
        let par = T::from_f64(100.0);
        let dt = spec.step_years_dt;
        let sqrt_dt = dt.sqrt();
        let cpn = par * spec.coupon_rate_c * dt;

        // Period-1 lognormal short-rate states r_{1,u} and r_{1,d}
        let r1_u = spec.initial_short_rate_r0
            * (spec.short_rate_drift_mu * dt
                + spec.lognormal_short_rate_volatility_sigma * sqrt_dt)
                .exp();
        let r1_d = spec.initial_short_rate_r0
            * (spec.short_rate_drift_mu * dt
                - spec.lognormal_short_rate_volatility_sigma * sqrt_dt)
                .exp();

        let disc0 = T::one()
            / (T::one() + (spec.initial_short_rate_r0 + spec.option_adjusted_spread_oas) * dt);
        let disc1_u = T::one() / (T::one() + (r1_u + spec.option_adjusted_spread_oas) * dt);
        let disc1_d = T::one() / (T::one() + (r1_d + spec.option_adjusted_spread_oas) * dt);

        // Straight ex-coupon continuation values at t = dt
        let v1_u_str = (par + cpn) * disc1_u;
        let v1_d_str = (par + cpn) * disc1_d;

        let half = T::from_f64(0.5);
        let v0_str = disc0 * (cpn + half * (v1_u_str + v1_d_str));

        // Smooth call exercise payoff max(V_1^str - K_call, 0) at t = dt
        let eps_sm = T::from_f64(1e-6);
        let diff_u = v1_u_str - spec.call_strike_price_k;
        let diff_d = v1_d_str - spec.call_strike_price_k;
        let call_intr_u = half * (diff_u + (diff_u * diff_u + eps_sm).sqrt());
        let call_intr_d = half * (diff_d + (diff_d * diff_d + eps_sm).sqrt());

        let v1_u_call = v1_u_str - call_intr_u;
        let v1_d_call = v1_d_str - call_intr_d;
        let v0_call = disc0 * (cpn + half * (v1_u_call + v1_d_call));

        let v_opt = v0_str - v0_call;

        // Net refunding savings after issuance cost c_ref in states where call is exercised
        let net_diff_u = v1_u_str - spec.call_strike_price_k - spec.refunding_transaction_cost;
        let net_diff_d = v1_d_str - spec.call_strike_price_k - spec.refunding_transaction_cost;
        let net_sav_u = half * (net_diff_u + (net_diff_u * net_diff_u + eps_sm).sqrt());
        let net_sav_d = half * (net_diff_d + (net_diff_d * net_diff_d + eps_sm).sqrt());
        let pv_net_savings = disc0 * half * (net_sav_u + net_sav_d);

        let efficiency = pv_net_savings / (v_opt + T::from_f64(1e-8));

        KalotayWilliamsFabozziCallableAnalytics {
            straight_option_free_bond_value: v0_str,
            callable_bond_value: v0_call,
            embedded_call_option_value: v_opt,
            refunding_efficiency_ratio: efficiency,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn kalotay_williams_fabozzi_callable_invariants_and_ad_parity() {
            let base = KalotayWilliamsFabozziCallableSpec {
                initial_short_rate_r0: 0.040_f64,
                short_rate_drift_mu: 0.002_f64,
                lognormal_short_rate_volatility_sigma: 0.22_f64,
                option_adjusted_spread_oas: 0.010_f64,
                coupon_rate_c: 0.058_f64,
                call_strike_price_k: 100.2_f64,
                refunding_transaction_cost: 0.35_f64,
                step_years_dt: 1.0_f64,
            };
            let res = evaluate_kalotay_williams_fabozzi_callable_bond_refunding_efficiency(&base);
            assert!(res.straight_option_free_bond_value > res.callable_bond_value);
            assert!(res.embedded_call_option_value > 0.0);
            assert!(res.refunding_efficiency_ratio > 0.0 && res.refunding_efficiency_ratio < 1.0);

            let inputs = [
                base.initial_short_rate_r0,
                base.short_rate_drift_mu,
                base.lognormal_short_rate_volatility_sigma,
                base.option_adjusted_spread_oas,
                base.coupon_rate_c,
                base.call_strike_price_k,
                base.refunding_transaction_cost,
                base.step_years_dt,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = KalotayWilliamsFabozziCallableSpec {
                initial_short_rate_r0: dual_in[0],
                short_rate_drift_mu: dual_in[1],
                lognormal_short_rate_volatility_sigma: dual_in[2],
                option_adjusted_spread_oas: dual_in[3],
                coupon_rate_c: dual_in[4],
                call_strike_price_k: dual_in[5],
                refunding_transaction_cost: dual_in[6],
                step_years_dt: dual_in[7],
            };
            let r_ad =
                evaluate_kalotay_williams_fabozzi_callable_bond_refunding_efficiency(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| KalotayWilliamsFabozziCallableSpec {
                    initial_short_rate_r0: v[0],
                    short_rate_drift_mu: v[1],
                    lognormal_short_rate_volatility_sigma: v[2],
                    option_adjusted_spread_oas: v[3],
                    coupon_rate_c: v[4],
                    call_strike_price_k: v[5],
                    refunding_transaction_cost: v[6],
                    step_years_dt: v[7],
                };
                let rp =
                    evaluate_kalotay_williams_fabozzi_callable_bond_refunding_efficiency(&mk(p));
                let rm =
                    evaluate_kalotay_williams_fabozzi_callable_bond_refunding_efficiency(&mk(m));
                let fd =
                    (rp.refunding_efficiency_ratio - rm.refunding_efficiency_ratio) / (2.0 * eps);
                assert!(
                    r_ad.refunding_efficiency_ratio.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.refunding_efficiency_ratio.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use kalotay_williams_fabozzi_callable_mod::*;
mod fabozzi_effective_duration_oas_mod {
    use ql_ad::Real;

    /// Input specification for the Fabozzi (2012, *The Handbook of Fixed Income Securities*,
    /// 8th ed., Ch. 38–40) & Tuckman & Serrat (2022, Ch. 19) Option-Adjusted Spread (OAS),
    /// effective duration $D_{\text{eff}}$, and negative effective convexity $\text{Conv}_{\text{eff}}$ engine.
    #[derive(Debug, Clone, Copy)]
    pub struct FabozziEffectiveDurationOasSpec<T: Real> {
        /// Benchmark risk-free yield $y_0$ (e.g., $0.042$).
        pub benchmark_yield_y0: T,
        /// Option-Adjusted Spread $\text{OAS}$ over the benchmark curve (e.g., $0.0085$).
        pub option_adjusted_spread_oas: T,
        /// Annual bond coupon rate $c$ per unit par (e.g., $0.056$).
        pub coupon_rate_c: T,
        /// Issuer call price $K_{\text{call}}$ per $100$ par (e.g., $100.5$).
        pub call_price_k: T,
        /// Forward price volatility $\sigma_P$ of the underlying option-free bond at call horizon (e.g., $0.065$).
        pub bond_price_volatility_sigma: T,
        /// Parallel curve shift $\Delta y > 0$ used for effective duration/convexity (e.g., $0.0025$).
        pub parallel_yield_shift_dy: T,
        /// Effective bond maturity $T$ in years (e.g., $5.0$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_fabozzi_tuckman_effective_duration_and_convexity_oas`].
    #[derive(Debug, Clone, Copy)]
    pub struct FabozziEffectiveDurationOasAnalytics<T: Real> {
        /// Option-free corporate bond value $V_{\text{nc}}(y_0 + \text{OAS})$ per $100$ par.
        pub option_free_bond_price: T,
        /// Callable corporate bond value $V_{\text{call}}(y_0 + \text{OAS}) = V_{\text{nc}} - C_{\text{call}}$.
        pub callable_bond_price: T,
        /// Fabozzi–Tuckman effective duration $D_{\text{eff}} = \frac{V_{\text{call}}(-\Delta y) - V_{\text{call}}(+\Delta y)}{2 V_{\text{call}}(0)\Delta y}$.
        pub effective_duration: T,
        /// Fabozzi–Tuckman effective convexity $\text{Conv}_{\text{eff}} = \frac{V_{\text{call}}(-\Delta y) + V_{\text{call}}(+\Delta y) - 2 V_{\text{call}}(0)}{V_{\text{call}}(0)(\Delta y)^2}$.
        pub effective_convexity: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    #[inline]
    fn std_norm_pdf<T: Real>(x: T) -> T {
        T::from_f64(0.398_942_280_401_432_7) * (-T::from_f64(0.5) * x * x).exp()
    }

    #[inline]
    fn callable_price_at_yield<T: Real>(y_tot: T, c: T, k_call: T, sig_p: T, t: T) -> (T, T) {
        let par = T::from_f64(100.0);
        let df = (-y_tot * t).exp();
        let annuity = (T::one() - df) / y_tot;
        let v_nc = par * (c * annuity + df);

        // Embedded issuer call option value via Bachelier normal bond-price call at half-maturity
        let t_call = T::from_f64(0.5) * t;
        let df_call = (-y_tot * t_call).exp();
        let std_p = par * sig_p * t_call.sqrt();
        let diff = v_nc - k_call;
        let z = diff / std_p;
        let call_val = df_call * (diff * std_norm_cdf(z) + std_p * std_norm_pdf(z));
        (v_nc, v_nc - call_val)
    }

    /// Evaluates the Fabozzi (2012, Ch. 38–40) & Tuckman & Serrat (2022, Ch. 19) callable bond
    /// price, effective duration $D_{\text{eff}}$, and negative effective convexity $\text{Conv}_{\text{eff}}$
    /// holding the Option-Adjusted Spread (OAS) constant across $\pm \Delta y$ curve shifts.
    pub fn evaluate_fabozzi_tuckman_effective_duration_and_convexity_oas<T: Real>(
        spec: &FabozziEffectiveDurationOasSpec<T>,
    ) -> FabozziEffectiveDurationOasAnalytics<T> {
        let y_base = spec.benchmark_yield_y0 + spec.option_adjusted_spread_oas;
        let dy = spec.parallel_yield_shift_dy;
        let c = spec.coupon_rate_c;
        let k_call = spec.call_price_k;
        let sig_p = spec.bond_price_volatility_sigma;
        let t = spec.maturity_years_t;

        let (v_nc_0, v_call_0) = callable_price_at_yield(y_base, c, k_call, sig_p, t);
        let (_, v_call_dn) = callable_price_at_yield(y_base - dy, c, k_call, sig_p, t);
        let (_, v_call_up) = callable_price_at_yield(y_base + dy, c, k_call, sig_p, t);

        let eff_dur = (v_call_dn - v_call_up) / (T::from_f64(2.0) * v_call_0 * dy);
        let eff_conv = (v_call_dn + v_call_up - T::from_f64(2.0) * v_call_0) / (v_call_0 * dy * dy);

        FabozziEffectiveDurationOasAnalytics {
            option_free_bond_price: v_nc_0,
            callable_bond_price: v_call_0,
            effective_duration: eff_dur,
            effective_convexity: eff_conv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn fabozzi_effective_duration_oas_invariants_and_ad_parity() {
            let base = FabozziEffectiveDurationOasSpec {
                benchmark_yield_y0: 0.042_f64,
                option_adjusted_spread_oas: 0.0085_f64,
                coupon_rate_c: 0.056_f64,
                call_price_k: 100.5_f64,
                bond_price_volatility_sigma: 0.065_f64,
                parallel_yield_shift_dy: 0.0025_f64,
                maturity_years_t: 5.0_f64,
            };
            let res = evaluate_fabozzi_tuckman_effective_duration_and_convexity_oas(&base);
            assert!(res.option_free_bond_price > res.callable_bond_price);
            assert!(res.effective_duration > 0.0 && res.effective_duration < base.maturity_years_t);
            // Near-the-money callable bond exhibits negative effective convexity
            assert!(res.effective_convexity < 0.0);

            let inputs = [
                base.benchmark_yield_y0,
                base.option_adjusted_spread_oas,
                base.coupon_rate_c,
                base.call_price_k,
                base.bond_price_volatility_sigma,
                base.parallel_yield_shift_dy,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = FabozziEffectiveDurationOasSpec {
                benchmark_yield_y0: dual_in[0],
                option_adjusted_spread_oas: dual_in[1],
                coupon_rate_c: dual_in[2],
                call_price_k: dual_in[3],
                bond_price_volatility_sigma: dual_in[4],
                parallel_yield_shift_dy: dual_in[5],
                maturity_years_t: dual_in[6],
            };
            let r_ad = evaluate_fabozzi_tuckman_effective_duration_and_convexity_oas(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| FabozziEffectiveDurationOasSpec {
                    benchmark_yield_y0: v[0],
                    option_adjusted_spread_oas: v[1],
                    coupon_rate_c: v[2],
                    call_price_k: v[3],
                    bond_price_volatility_sigma: v[4],
                    parallel_yield_shift_dy: v[5],
                    maturity_years_t: v[6],
                };
                let rp = evaluate_fabozzi_tuckman_effective_duration_and_convexity_oas(&mk(p));
                let rm = evaluate_fabozzi_tuckman_effective_duration_and_convexity_oas(&mk(m));
                let fd = (rp.effective_duration - rm.effective_duration) / (2.0 * eps);
                assert!(
                    r_ad.effective_duration.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.effective_duration.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use fabozzi_effective_duration_oas_mod::*;
mod duffie_krishnamurthy_special_repo_mod {
    use ql_ad::Real;

    /// Input specification for the Duffie (1996, *Journal of Finance* 51(2), *"Special Repo Rates"*)
    /// & Krishnamurthy (2002, *Journal of Financial Economics* 66(2)) On-the-Run vs. Off-the-Run
    /// Treasury special repo financing dividend and auction-cycle convergence engine.
    #[derive(Debug, Clone, Copy)]
    pub struct DuffieKrishnamurthySpecialRepoTreasurySpec<T: Real> {
        /// General Collateral (GC) repo rate $r_{\text{GC}}$ (e.g., $0.045$).
        pub general_collateral_repo_rate: T,
        /// Initial on-the-run repo specialness spread $s_0 = r_{\text{GC}} - r_{\text{spec}}(0)$ (e.g., $0.0110$).
        pub initial_repo_specialness_s0: T,
        /// Auction-cycle specialness decay rate $\kappa_s > 0$ (e.g., $3.20$).
        pub specialness_decay_rate_kappa: T,
        /// Short-interest squeeze sensitivity $\gamma_{\text{short}} > 0$ (e.g., $0.0085$).
        pub short_squeeze_sensitivity_gamma: T,
        /// Short-interest to floating supply ratio $q_{\text{short}} \in (0, 1)$ (e.g., $0.35$).
        pub short_to_float_ratio_q: T,
        /// Time to next benchmark Treasury auction $T_{\text{auc}}$ in years (e.g., $0.25$).
        pub remaining_auction_cycle_years: T,
        /// Benchmark Treasury modified duration $D_{\text{mod}}$ in years (e.g., $7.80$).
        pub modified_duration_d_mod: T,
        /// Dealer repo intermediation haircut fraction $h_{\text{dealer}} \in (0, 1)$ (e.g., $0.18$).
        pub dealer_intermediation_haircut_h: T,
    }

    /// Analytics returned by [`evaluate_duffie_krishnamurthy_on_the_run_special_repo_premium`].
    #[derive(Debug, Clone, Copy)]
    pub struct DuffieKrishnamurthySpecialRepoTreasuryAnalytics<T: Real> {
        /// Present value of special repo financing dividends $\Delta P_{\text{OTR}}$ per $100$ par.
        pub on_the_run_special_repo_price_premium: T,
        /// Equilibrium On-the-Run vs. Off-the-Run yield spread compression $\Delta y_{\text{OTR-OFR}}$ (in bps).
        pub on_the_run_yield_spread_compression_bps: T,
        /// Net convergence arbitrage carry per $100$ par after dealer repo intermediation haircut.
        pub net_convergence_arbitrage_carry: T,
        /// Effective annualized duration-hedged specialness yield pick-up (in bps/yr).
        pub annualized_specialness_pickup_rate: T,
    }

    /// Evaluates the Duffie (1996, Proposition 1) & Krishnamurthy (2002, Eqs. 3–8) present value
    /// of special repo financing cashflows over the auction cycle $[0, T_{\text{auc}}]$:
    /// $$\Delta P_{\text{OTR}} = 100 \left( s_0 \frac{1 - e^{-(r_{\text{GC}} + \kappa_s)T_{\text{auc}}}}{r_{\text{GC}} + \kappa_s} + \gamma_{\text{short}} q_{\text{short}} \frac{1 - e^{-(r_{\text{GC}} + \tfrac{1}{2}\kappa_s)T_{\text{auc}}}}{r_{\text{GC}} + \tfrac{1}{2}\kappa_s} \right),$$
    /// and the implied yield spread compression $\Delta y_{\text{OTR-OFR}}^{\text{bps}} = 100 \cdot \Delta P_{\text{OTR}} / D_{\text{mod}}$.
    pub fn evaluate_duffie_krishnamurthy_on_the_run_special_repo_premium<T: Real>(
        spec: &DuffieKrishnamurthySpecialRepoTreasurySpec<T>,
    ) -> DuffieKrishnamurthySpecialRepoTreasuryAnalytics<T> {
        let r_gc = spec.general_collateral_repo_rate;
        let s0 = spec.initial_repo_specialness_s0;
        let kap = spec.specialness_decay_rate_kappa;
        let gam = spec.short_squeeze_sensitivity_gamma;
        let q_sh = spec.short_to_float_ratio_q;
        let t_auc = spec.remaining_auction_cycle_years;
        let d_mod = spec.modified_duration_d_mod;
        let h_dlr = spec.dealer_intermediation_haircut_h;

        let rate1 = r_gc + kap;
        let rate2 = r_gc + T::from_f64(0.5) * kap;
        let pv1 = s0 * (T::one() - (-rate1 * t_auc).exp()) / rate1;
        let pv2 = gam * q_sh * (T::one() - (-rate2 * t_auc).exp()) / rate2;

        let dp_otr = T::from_f64(100.0) * (pv1 + pv2);
        let dy_bps = (dp_otr / d_mod) * T::from_f64(100.0);
        let net_carry = (T::one() - h_dlr) * dp_otr;
        let ann_pickup = (T::one() - h_dlr) * dy_bps / t_auc;

        DuffieKrishnamurthySpecialRepoTreasuryAnalytics {
            on_the_run_special_repo_price_premium: dp_otr,
            on_the_run_yield_spread_compression_bps: dy_bps,
            net_convergence_arbitrage_carry: net_carry,
            annualized_specialness_pickup_rate: ann_pickup,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn duffie_krishnamurthy_special_repo_invariants_and_ad_parity() {
            let base = DuffieKrishnamurthySpecialRepoTreasurySpec {
                general_collateral_repo_rate: 0.045_f64,
                initial_repo_specialness_s0: 0.0110_f64,
                specialness_decay_rate_kappa: 3.20_f64,
                short_squeeze_sensitivity_gamma: 0.0085_f64,
                short_to_float_ratio_q: 0.35_f64,
                remaining_auction_cycle_years: 0.25_f64,
                modified_duration_d_mod: 7.80_f64,
                dealer_intermediation_haircut_h: 0.18_f64,
            };
            let res = evaluate_duffie_krishnamurthy_on_the_run_special_repo_premium(&base);
            assert!(res.on_the_run_special_repo_price_premium > 0.0);
            assert!(res.on_the_run_yield_spread_compression_bps > 0.0);
            assert!(
                res.net_convergence_arbitrage_carry < res.on_the_run_special_repo_price_premium
                    && res.net_convergence_arbitrage_carry > 0.0
            );
            assert!(res.annualized_specialness_pickup_rate > 0.0);

            let inputs = [
                base.general_collateral_repo_rate,
                base.initial_repo_specialness_s0,
                base.specialness_decay_rate_kappa,
                base.short_squeeze_sensitivity_gamma,
                base.short_to_float_ratio_q,
                base.remaining_auction_cycle_years,
                base.modified_duration_d_mod,
                base.dealer_intermediation_haircut_h,
            ];
            let mut dual_in = [Dual::<8>::from_value(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = DuffieKrishnamurthySpecialRepoTreasurySpec {
                general_collateral_repo_rate: dual_in[0],
                initial_repo_specialness_s0: dual_in[1],
                specialness_decay_rate_kappa: dual_in[2],
                short_squeeze_sensitivity_gamma: dual_in[3],
                short_to_float_ratio_q: dual_in[4],
                remaining_auction_cycle_years: dual_in[5],
                modified_duration_d_mod: dual_in[6],
                dealer_intermediation_haircut_h: dual_in[7],
            };
            let r_ad = evaluate_duffie_krishnamurthy_on_the_run_special_repo_premium(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| DuffieKrishnamurthySpecialRepoTreasurySpec {
                    general_collateral_repo_rate: v[0],
                    initial_repo_specialness_s0: v[1],
                    specialness_decay_rate_kappa: v[2],
                    short_squeeze_sensitivity_gamma: v[3],
                    short_to_float_ratio_q: v[4],
                    remaining_auction_cycle_years: v[5],
                    modified_duration_d_mod: v[6],
                    dealer_intermediation_haircut_h: v[7],
                };
                let rp = evaluate_duffie_krishnamurthy_on_the_run_special_repo_premium(&mk(p));
                let rm = evaluate_duffie_krishnamurthy_on_the_run_special_repo_premium(&mk(m));
                let fd = (rp.annualized_specialness_pickup_rate
                    - rm.annualized_specialness_pickup_rate)
                    / (2.0 * eps);
                assert!(
                    r_ad.annualized_specialness_pickup_rate.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.annualized_specialness_pickup_rate.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use duffie_krishnamurthy_special_repo_mod::*;

/// Input specification for the Acharya, Amihud & Bharath (2013, *Journal of
/// Financial Economics* 110(2)) Stress-Regime Flight-to-Liquidity Corporate Bond
/// Beta & Spread engine.
#[derive(Debug, Clone, Copy)]
pub struct AcharyaAmihudBharathLiquidityStressBondSpec<T: Real> {
    /// Baseline structural default credit spread $s_{\text{def}} > 0$.
    pub structural_default_spread_s_def: T,
    /// Amihud/bid-ask bond illiquidity level $\text{ILLIQ}_i > 0$.
    pub bond_illiquidity_level_illiq: T,
    /// Normal-regime liquidity beta $\beta_{i,\text{liq}}^{(1)}$.
    pub normal_regime_liquidity_beta_1: T,
    /// Incremental stress-regime liquidity beta $\beta_{i,\text{liq}}^{(2)}$.
    pub stress_regime_liquidity_beta_2: T,
    /// Filtered Markov probability of financial stress regime $\pi_S \in (0,1)$.
    pub stress_regime_probability_pi_s: T,
    /// Market price of aggregate illiquidity risk $\lambda_{\text{liq}} > 0$.
    pub market_illiquidity_risk_price_lambda_l: T,
    /// Flight-to-quality spread compression wedge $\psi_{\text{FTQ}}$.
    pub flight_to_quality_treasury_wedge_psi: T,
    /// Corporate bond modified duration $D_{\text{mod}} > 0$.
    pub modified_duration_d: T,
}

/// Analytics returned by [`evaluate_acharya_amihud_bharath_liquidity_stress_bond_spread`].
#[derive(Debug, Clone, Copy)]
pub struct AcharyaAmihudBharathLiquidityStressBondAnalytics<T: Real> {
    /// Effective regime-weighted liquidity beta $\beta_{i,\text{liq}}^{\text{eff}} = \beta_{i,\text{liq}}^{(1)} + \pi_S\beta_{i,\text{liq}}^{(2)}$.
    pub regime_weighted_liquidity_beta: T,
    /// Liquidity risk premium spread component $s_{\text{liq}}$.
    pub liquidity_risk_premium_spread: T,
    /// Total stress-adjusted corporate bond yield spread $s_{\text{total}} = s_{\text{def}} + s_{\text{liq}}$.
    pub stress_adjusted_corporate_bond_spread: T,
    /// Expected bond excess return shock under stress regime transition $\Delta R_{\text{stress}}$.
    pub expected_stress_return_impact: T,
}

/// Evaluates the Acharya, Amihud & Bharath (2013, *JFE* 110(2)) regime-switching
/// corporate bond liquidity beta, flight-to-quality spread bifurcation, and
/// stress-regime corporate yield spread.
pub fn evaluate_acharya_amihud_bharath_liquidity_stress_bond_spread<T: Real>(
    spec: &AcharyaAmihudBharathLiquidityStressBondSpec<T>,
) -> AcharyaAmihudBharathLiquidityStressBondAnalytics<T> {
    let s_def = spec
        .structural_default_spread_s_def
        .max_val(T::from_f64(1e-5));
    let illiq = spec.bond_illiquidity_level_illiq.max_val(T::from_f64(1e-5));
    let b1 = spec.normal_regime_liquidity_beta_1;
    let b2 = spec.stress_regime_liquidity_beta_2;
    let pi_s = spec
        .stress_regime_probability_pi_s
        .min_val(T::from_f64(0.999))
        .max_val(T::from_f64(1e-4));
    let lam_l = spec
        .market_illiquidity_risk_price_lambda_l
        .max_val(T::from_f64(1e-5));
    let psi = spec.flight_to_quality_treasury_wedge_psi;
    let d_mod = spec.modified_duration_d.max_val(T::from_f64(0.2));

    let beta_eff = b1 + pi_s * b2;
    let ftq_term = psi * (-illiq).exp();
    let s_liq = (illiq / d_mod) * (T::one() + T::from_f64(5.0) * s_def) + lam_l * beta_eff
        - pi_s * ftq_term;
    let s_tot = s_def + s_liq;
    let stress_ret = -d_mod * (lam_l * b2 - ftq_term + illiq * s_def) * pi_s;

    AcharyaAmihudBharathLiquidityStressBondAnalytics {
        regime_weighted_liquidity_beta: beta_eff,
        liquidity_risk_premium_spread: s_liq,
        stress_adjusted_corporate_bond_spread: s_tot,
        expected_stress_return_impact: stress_ret,
    }
}

#[cfg(test)]
mod acharya_amihud_bharath_liquidity_stress_bond_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn acharya_amihud_bharath_liquidity_stress_bond_invariants_and_ad_parity() {
        let base = AcharyaAmihudBharathLiquidityStressBondSpec {
            structural_default_spread_s_def: 0.018_f64,
            bond_illiquidity_level_illiq: 0.012_f64,
            normal_regime_liquidity_beta_1: 0.25_f64,
            stress_regime_liquidity_beta_2: 0.85_f64,
            stress_regime_probability_pi_s: 0.30_f64,
            market_illiquidity_risk_price_lambda_l: 0.015_f64,
            flight_to_quality_treasury_wedge_psi: 0.004_f64,
            modified_duration_d: 6.2_f64,
        };
        let r0 = evaluate_acharya_amihud_bharath_liquidity_stress_bond_spread(&base);
        assert!(r0.regime_weighted_liquidity_beta > base.normal_regime_liquidity_beta_1);
        assert!(r0.stress_adjusted_corporate_bond_spread > base.structural_default_spread_s_def);
        assert!(r0.expected_stress_return_impact < 0.0);

        let x0 = [
            base.structural_default_spread_s_def,
            base.bond_illiquidity_level_illiq,
            base.normal_regime_liquidity_beta_1,
            base.stress_regime_liquidity_beta_2,
            base.stress_regime_probability_pi_s,
            base.market_illiquidity_risk_price_lambda_l,
            base.flight_to_quality_treasury_wedge_psi,
            base.modified_duration_d,
        ];
        let mut xd = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_acharya_amihud_bharath_liquidity_stress_bond_spread(
            &AcharyaAmihudBharathLiquidityStressBondSpec {
                structural_default_spread_s_def: xd[0],
                bond_illiquidity_level_illiq: xd[1],
                normal_regime_liquidity_beta_1: xd[2],
                stress_regime_liquidity_beta_2: xd[3],
                stress_regime_probability_pi_s: xd[4],
                market_illiquidity_risk_price_lambda_l: xd[5],
                flight_to_quality_treasury_wedge_psi: xd[6],
                modified_duration_d: xd[7],
            },
        );
        let h = 1e-6;
        for i in 0..8 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_acharya_amihud_bharath_liquidity_stress_bond_spread(
                &AcharyaAmihudBharathLiquidityStressBondSpec {
                    structural_default_spread_s_def: xp[0],
                    bond_illiquidity_level_illiq: xp[1],
                    normal_regime_liquidity_beta_1: xp[2],
                    stress_regime_liquidity_beta_2: xp[3],
                    stress_regime_probability_pi_s: xp[4],
                    market_illiquidity_risk_price_lambda_l: xp[5],
                    flight_to_quality_treasury_wedge_psi: xp[6],
                    modified_duration_d: xp[7],
                },
            );
            let rm = evaluate_acharya_amihud_bharath_liquidity_stress_bond_spread(
                &AcharyaAmihudBharathLiquidityStressBondSpec {
                    structural_default_spread_s_def: xm[0],
                    bond_illiquidity_level_illiq: xm[1],
                    normal_regime_liquidity_beta_1: xm[2],
                    stress_regime_liquidity_beta_2: xm[3],
                    stress_regime_probability_pi_s: xm[4],
                    market_illiquidity_risk_price_lambda_l: xm[5],
                    flight_to_quality_treasury_wedge_psi: xm[6],
                    modified_duration_d: xm[7],
                },
            );
            let fd = (rp.stress_adjusted_corporate_bond_spread
                - rm.stress_adjusted_corporate_bond_spread)
                / (2.0 * h);
            assert!(rad.stress_adjusted_corporate_bond_spread.d[i].abs() > 1e-9);
            assert!((rad.stress_adjusted_corporate_bond_spread.d[i] - fd).abs() < 1e-6);
        }
    }
}

/// Input specification for the Longstaff, Mithal & Neis (2005, *Journal of
/// Finance* 60(5), "Corporate Yield Spreads: Default Risk or Liquidity? New
/// Evidence from the Credit Default Swap Market") CDS-Bond Basis Non-Default
/// Liquidity & Tax-Coupon Wedge engine.
#[derive(Debug, Clone, Copy)]
pub struct LongstaffMithalNeisCdsBondLiquidityBasisSpec<T: Real> {
    /// Observed corporate bond yield $y_{\text{corp}} > 0$.
    pub corporate_bond_yield_y_corp: T,
    /// Maturity-matched Treasury benchmark yield $y_{\text{tsy}} > 0$.
    pub treasury_benchmark_yield_y_tsy: T,
    /// Pure default-risk CDS spread $s_{\text{CDS}} > 0$.
    pub cds_default_spread_s_cds: T,
    /// Corporate bond coupon rate $c_{\text{corp}} > 0$.
    pub corporate_coupon_rate_c: T,
    /// Marginal state/local tax differential on corporate vs. Treasury coupons $\tau_{\text{tax}} \in (0, 0.20)$.
    pub state_tax_differential_tau: T,
    /// Corporate bond bid-ask illiquidity spread $s_{\text{bid-ask}} > 0$.
    pub bond_bid_ask_spread_s_ba: T,
    /// Outstanding bond issue size / turnover liquidity factor $\ell_{\text{size}} > 0$.
    pub issue_turnover_liquidity_factor_l: T,
    /// Bond modified duration $D_{\text{mod}} > 0$.
    pub modified_duration_d: T,
}

/// Analytics returned by [`evaluate_longstaff_mithal_neis_cds_bond_liquidity_basis`].
#[derive(Debug, Clone, Copy)]
pub struct LongstaffMithalNeisCdsBondLiquidityBasisAnalytics<T: Real> {
    /// Total observed corporate yield spread over Treasuries $S_{\text{total}} = y_{\text{corp}} - y_{\text{tsy}}$.
    pub total_corporate_treasury_spread: T,
    /// Default-related fraction of corporate spread $s_{\text{CDS}} / S_{\text{total}}$.
    pub default_component_share: T,
    /// Non-default illiquidity and tax-asymmetry yield wedge $S_{\text{non-def}}$ in basis points.
    pub non_default_liquidity_tax_wedge_bps: T,
    /// Implied negative CDS-bond basis $(s_{\text{CDS}} - S_{\text{corp,adj}})$ in basis points.
    pub adjusted_cds_bond_basis_bps: T,
}

/// Evaluates the Longstaff, Mithal & Neis (2005, *JF* 60(5)) decomposition of
/// corporate bond yield spreads into pure default risk (spanned by CDS spreads),
/// state/local coupon tax asymmetry, and bond-specific market illiquidity.
pub fn evaluate_longstaff_mithal_neis_cds_bond_liquidity_basis<T: Real>(
    spec: &LongstaffMithalNeisCdsBondLiquidityBasisSpec<T>,
) -> LongstaffMithalNeisCdsBondLiquidityBasisAnalytics<T> {
    let one = T::one();
    let y_corp = spec.corporate_bond_yield_y_corp;
    let y_tsy = spec.treasury_benchmark_yield_y_tsy;
    let s_cds = spec.cds_default_spread_s_cds.max_val(T::from_f64(1e-5));
    let c_corp = spec.corporate_coupon_rate_c.max_val(T::from_f64(1e-4));
    let tau_tax = spec
        .state_tax_differential_tau
        .min_val(T::from_f64(0.25))
        .max_val(T::from_f64(1e-4));
    let s_ba = spec.bond_bid_ask_spread_s_ba.max_val(T::from_f64(1e-5));
    let l_turn = spec
        .issue_turnover_liquidity_factor_l
        .max_val(T::from_f64(0.05));
    let d_mod = spec.modified_duration_d.max_val(T::from_f64(0.5));

    let s_total = (y_corp - y_tsy).max_val(T::from_f64(1e-4));
    let def_share = (s_cds / s_total).min_val(one);

    // Tax wedge from state/local taxation of corporate coupons relative to exempt Treasuries
    let tax_wedge = tau_tax * c_corp;
    // Illiquidity yield component amortized over duration and dampened by turnover
    let liq_wedge = (s_ba / d_mod) * (one + T::from_f64(0.5) / l_turn);
    let non_def_bps = (tax_wedge + liq_wedge) * T::from_f64(10_000.0);

    let basis_bps = (s_cds - (s_total - tax_wedge - liq_wedge)) * T::from_f64(10_000.0);

    LongstaffMithalNeisCdsBondLiquidityBasisAnalytics {
        total_corporate_treasury_spread: s_total,
        default_component_share: def_share,
        non_default_liquidity_tax_wedge_bps: non_def_bps,
        adjusted_cds_bond_basis_bps: basis_bps,
    }
}

#[cfg(test)]
mod longstaff_mithal_neis_cds_bond_basis_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn longstaff_mithal_neis_cds_bond_basis_invariants_and_ad_parity() {
        let base = LongstaffMithalNeisCdsBondLiquidityBasisSpec {
            corporate_bond_yield_y_corp: 0.058_f64,
            treasury_benchmark_yield_y_tsy: 0.040_f64,
            cds_default_spread_s_cds: 0.0115_f64,
            corporate_coupon_rate_c: 0.055_f64,
            state_tax_differential_tau: 0.065_f64,
            bond_bid_ask_spread_s_ba: 0.0048_f64,
            issue_turnover_liquidity_factor_l: 1.25_f64,
            modified_duration_d: 6.0_f64,
        };
        let r0 = evaluate_longstaff_mithal_neis_cds_bond_liquidity_basis(&base);
        // LMN (2005): default component is the majority (>50%) of corporate spread, with positive non-default liquidity/tax wedge
        assert!(r0.default_component_share > 0.50 && r0.default_component_share < 1.0);
        assert!(r0.non_default_liquidity_tax_wedge_bps > 0.0);

        let x0 = [
            base.corporate_bond_yield_y_corp,
            base.treasury_benchmark_yield_y_tsy,
            base.cds_default_spread_s_cds,
            base.corporate_coupon_rate_c,
            base.state_tax_differential_tau,
            base.bond_bid_ask_spread_s_ba,
            base.issue_turnover_liquidity_factor_l,
            base.modified_duration_d,
        ];
        let mut xd = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_longstaff_mithal_neis_cds_bond_liquidity_basis(
            &LongstaffMithalNeisCdsBondLiquidityBasisSpec {
                corporate_bond_yield_y_corp: xd[0],
                treasury_benchmark_yield_y_tsy: xd[1],
                cds_default_spread_s_cds: xd[2],
                corporate_coupon_rate_c: xd[3],
                state_tax_differential_tau: xd[4],
                bond_bid_ask_spread_s_ba: xd[5],
                issue_turnover_liquidity_factor_l: xd[6],
                modified_duration_d: xd[7],
            },
        );
        let h = 1e-6;
        for i in 0..8 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_longstaff_mithal_neis_cds_bond_liquidity_basis(
                &LongstaffMithalNeisCdsBondLiquidityBasisSpec {
                    corporate_bond_yield_y_corp: xp[0],
                    treasury_benchmark_yield_y_tsy: xp[1],
                    cds_default_spread_s_cds: xp[2],
                    corporate_coupon_rate_c: xp[3],
                    state_tax_differential_tau: xp[4],
                    bond_bid_ask_spread_s_ba: xp[5],
                    issue_turnover_liquidity_factor_l: xp[6],
                    modified_duration_d: xp[7],
                },
            );
            let rm = evaluate_longstaff_mithal_neis_cds_bond_liquidity_basis(
                &LongstaffMithalNeisCdsBondLiquidityBasisSpec {
                    corporate_bond_yield_y_corp: xm[0],
                    treasury_benchmark_yield_y_tsy: xm[1],
                    cds_default_spread_s_cds: xm[2],
                    corporate_coupon_rate_c: xm[3],
                    state_tax_differential_tau: xm[4],
                    bond_bid_ask_spread_s_ba: xm[5],
                    issue_turnover_liquidity_factor_l: xm[6],
                    modified_duration_d: xm[7],
                },
            );
            let fd = (rp.adjusted_cds_bond_basis_bps - rm.adjusted_cds_bond_basis_bps) / (2.0 * h);
            assert!(rad.adjusted_cds_bond_basis_bps.d[i].abs() > 1e-9);
            assert!((rad.adjusted_cds_bond_basis_bps.d[i] - fd).abs() < 1e-4);
        }
    }
}

/// Input specification for the Bao, Pan & Wang (2011, *Journal of Finance* 66(3),
/// "The Illiquidity of Corporate Bonds") Roll-Style Negative Autocovariance
/// Corporate Bond Illiquidity Measure $\gamma_{\text{BPW}} = -\text{Cov}(\Delta p_t, \Delta p_{t-1})$
/// & Credit Yield Spread Wedge engine.
#[derive(Debug, Clone, Copy)]
pub struct BaoPanWangBondAutocovarianceIlliquiditySpec<T: Real> {
    /// Sample lag-1 log price change autocovariance $\text{Cov}(\Delta p_t, \Delta p_{t-1}) < 0$ from TRACE transaction data.
    pub lag1_log_price_autocovariance_cov1: T,
    /// Quoted dealer bid-ask spread $s_{\text{BA}} > 0$ (as a fraction of par).
    pub quoted_bid_ask_spread_s_ba: T,
    /// Aggregate market VIX / dealer funding stress index $\sigma_{\text{VIX}} > 0$.
    pub market_vix_stress_sigma_vix: T,
    /// Bond age / off-the-run seasoning in years $A_{\text{bond}} > 0$.
    pub bond_age_years_a: T,
    /// Baseline structural default yield spread $s_{\text{def}} > 0$.
    pub structural_default_spread_s_def: T,
    /// Cross-sectional yield-spread loading $\beta_\gamma > 0$ on $\gamma_{\text{BPW}}$.
    pub illiquidity_yield_loading_beta_gamma: T,
}

/// Analytics returned by [`evaluate_bao_pan_wang_corporate_bond_autocovariance_illiquidity`].
#[derive(Debug, Clone, Copy)]
pub struct BaoPanWangBondAutocovarianceIlliquidityAnalytics<T: Real> {
    /// Bao-Pan-Wang transient price-reversal illiquidity measure $\gamma_{\text{BPW}} = -\text{Cov}(\Delta p_t, \Delta p_{t-1})$.
    pub bao_pan_wang_gamma_illiquidity: T,
    /// Implied Roll effective round-trip transaction cost $2\sqrt{\gamma_{\text{BPW}}}$ (in bps).
    pub implied_effective_spread_bps: T,
    /// Ratio of $\gamma_{\text{BPW}}$ to quoted Roll bid-ask bounce $\frac{1}{4} s_{\text{BA}}^2$.
    pub excess_reversal_over_quoted_bid_ask_ratio: T,
    /// Total corporate bond yield spread $s_{\text{total}} = s_{\text{def}} + \Delta s_{\text{illiq}}$ (in bps).
    pub total_corporate_bond_spread_bps: T,
}

/// Evaluates the Bao, Pan & Wang (2011, *JF* 66(3)) corporate bond illiquidity
/// metric $\gamma_{\text{BPW}} = -\text{Cov}(\Delta p_t, \Delta p_{t-1})$ and
/// its amplification by VIX funding stress and bond seasoning:
/// $$\gamma_{\text{BPW}} = \left(-\text{Cov}(\Delta p_t, \Delta p_{t-1})\right) \left(1 + 0.8\,\sigma_{\text{VIX}} + 0.06\,A_{\text{bond}}\right),$$
/// quantifying the excess transient reversal over the quoted bid-ask benchmark
/// $\frac{1}{4}s_{\text{BA}}^2$ and the resulting corporate bond yield spread in bps.
pub fn evaluate_bao_pan_wang_corporate_bond_autocovariance_illiquidity<T: Real>(
    spec: &BaoPanWangBondAutocovarianceIlliquiditySpec<T>,
) -> BaoPanWangBondAutocovarianceIlliquidityAnalytics<T> {
    let one = T::one();
    let bps = T::from_f64(10_000.0);

    let neg_cov1 = (-spec.lag1_log_price_autocovariance_cov1).max_val(T::from_f64(1e-6));
    let s_ba = spec.quoted_bid_ask_spread_s_ba.max_val(T::from_f64(1e-4));
    let vix = spec.market_vix_stress_sigma_vix.max_val(T::from_f64(0.05));
    let age = spec.bond_age_years_a.max_val(T::from_f64(0.1));
    let s_def = spec
        .structural_default_spread_s_def
        .max_val(T::from_f64(1e-4));
    let beta_g = spec
        .illiquidity_yield_loading_beta_gamma
        .max_val(T::from_f64(0.1));

    let gamma_bpw = neg_cov1 * (one + T::from_f64(0.8) * vix + T::from_f64(0.06) * age);
    let eff_spread = T::from_f64(2.0) * gamma_bpw.sqrt();
    let roll_quoted_var = T::from_f64(0.25) * s_ba * s_ba;
    let excess_ratio = gamma_bpw / roll_quoted_var;

    let illiq_spread = beta_g * gamma_bpw + T::from_f64(0.35) * s_ba;
    let total_spread_bps = (s_def + illiq_spread) * bps;

    BaoPanWangBondAutocovarianceIlliquidityAnalytics {
        bao_pan_wang_gamma_illiquidity: gamma_bpw,
        implied_effective_spread_bps: eff_spread * bps,
        excess_reversal_over_quoted_bid_ask_ratio: excess_ratio,
        total_corporate_bond_spread_bps: total_spread_bps,
    }
}

#[cfg(test)]
mod bao_pan_wang_bond_illiquidity_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn bao_pan_wang_bond_illiquidity_invariants_and_ad_parity() {
        let base = BaoPanWangBondAutocovarianceIlliquiditySpec {
            lag1_log_price_autocovariance_cov1: -0.00011_f64,
            quoted_bid_ask_spread_s_ba: 0.0045_f64,
            market_vix_stress_sigma_vix: 0.22_f64,
            bond_age_years_a: 3.5_f64,
            structural_default_spread_s_def: 0.0095_f64,
            illiquidity_yield_loading_beta_gamma: 38.0_f64,
        };
        let r0 = evaluate_bao_pan_wang_corporate_bond_autocovariance_illiquidity(&base);
        // Bao, Pan & Wang (2011): empirical gamma_BPW substantially exceeds quoted Roll bid-ask bounce
        assert!(r0.bao_pan_wang_gamma_illiquidity > 0.0);
        assert!(r0.excess_reversal_over_quoted_bid_ask_ratio > 1.0);
        assert!(
            r0.total_corporate_bond_spread_bps > base.structural_default_spread_s_def * 10_000.0
        );

        let x0 = [
            base.lag1_log_price_autocovariance_cov1,
            base.quoted_bid_ask_spread_s_ba,
            base.market_vix_stress_sigma_vix,
            base.bond_age_years_a,
            base.structural_default_spread_s_def,
            base.illiquidity_yield_loading_beta_gamma,
        ];
        let mut xd = [Dual::<6>::from_f64(0.0); 6];
        for i in 0..6 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_bao_pan_wang_corporate_bond_autocovariance_illiquidity(
            &BaoPanWangBondAutocovarianceIlliquiditySpec {
                lag1_log_price_autocovariance_cov1: xd[0],
                quoted_bid_ask_spread_s_ba: xd[1],
                market_vix_stress_sigma_vix: xd[2],
                bond_age_years_a: xd[3],
                structural_default_spread_s_def: xd[4],
                illiquidity_yield_loading_beta_gamma: xd[5],
            },
        );
        let h = 1e-7;
        for i in 0..6 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_bao_pan_wang_corporate_bond_autocovariance_illiquidity(
                &BaoPanWangBondAutocovarianceIlliquiditySpec {
                    lag1_log_price_autocovariance_cov1: xp[0],
                    quoted_bid_ask_spread_s_ba: xp[1],
                    market_vix_stress_sigma_vix: xp[2],
                    bond_age_years_a: xp[3],
                    structural_default_spread_s_def: xp[4],
                    illiquidity_yield_loading_beta_gamma: xp[5],
                },
            );
            let rm = evaluate_bao_pan_wang_corporate_bond_autocovariance_illiquidity(
                &BaoPanWangBondAutocovarianceIlliquiditySpec {
                    lag1_log_price_autocovariance_cov1: xm[0],
                    quoted_bid_ask_spread_s_ba: xm[1],
                    market_vix_stress_sigma_vix: xm[2],
                    bond_age_years_a: xm[3],
                    structural_default_spread_s_def: xm[4],
                    illiquidity_yield_loading_beta_gamma: xm[5],
                },
            );
            let fd = (rp.total_corporate_bond_spread_bps - rm.total_corporate_bond_spread_bps)
                / (2.0 * h);
            assert!(rad.total_corporate_bond_spread_bps.d[i].abs() > 1e-9);
            assert!(
                (rad.total_corporate_bond_spread_bps.d[i] - fd).abs()
                    / rad.total_corporate_bond_spread_bps.d[i].abs().max(1.0)
                    < 1e-4
            );
        }
    }
}

/// Input specification for the Chen, Lesmond & Wei (2007, *Journal of Finance*
/// 62(1), 119–149, "Corporate Yield Spreads and Bond Liquidity") Lesmond–Ogden–
/// Trzcinka (LOT) Zero-Return Threshold Corporate Bond Illiquidity Model.
#[derive(Debug, Clone, Copy)]
pub struct ChenLesmondWeiLotBondIlliquiditySpec<T: Real> {
    /// Observed proportion of zero-return trading days $z_{\text{zero}} \in (0, 1)$.
    pub zero_return_day_proportion_z: T,
    /// Latent true daily bond return volatility $\sigma_{\text{latent}} > 0$.
    pub latent_daily_return_volatility_sigma: T,
    /// Quoted dealer round-trip bid-ask spread in bps $s_{\text{ba,bps}} > 0$.
    pub quoted_bid_ask_spread_bps: T,
    /// Bond maturity duration in years $D_{\text{mat}} > 0$.
    pub bond_maturity_years_d: T,
    /// Issuer debt-to-capital leverage ratio $L_{\text{lev}} \in (0, 1)$.
    pub issuer_leverage_ratio_l: T,
    /// Yield spread slope on LOT transaction costs $\beta_{\text{LOT}} > 0$.
    pub lot_yield_spread_loading_beta: T,
}

/// Analytics returned by [`evaluate_chen_lesmond_wei_lot_zero_return_bond_illiquidity`].
#[derive(Debug, Clone, Copy)]
pub struct ChenLesmondWeiLotBondIlliquidityAnalytics<T: Real> {
    /// Implied LOT round-trip transaction cost band $\text{LOT}_i = \alpha_{2,i} - \alpha_{1,i}$ (in bps).
    pub lot_roundtrip_transaction_cost_bps: T,
    /// Zero-return threshold width contribution $2\sigma_{\text{latent}}\Phi^{-1}((1+z_{\text{zero}})/2)$ (in bps).
    pub zero_return_threshold_band_bps: T,
    /// Liquidity component of the corporate bond yield spread (in bps).
    pub lot_corporate_yield_spread_component_bps: T,
}

/// Evaluates the Chen, Lesmond & Wei (2007, *JF* 62(1)) corporate bond
/// illiquidity model based on the Lesmond, Ogden & Trzcinka (1999) zero-return
/// threshold estimator $\text{LOT}_i = \alpha_{2,i} - \alpha_{1,i}$. Using the
/// inverse-normal Winitzki/probit expansion for the zero-return probability
/// $z_{\text{zero}} = \Phi(\alpha_{2,i}/\sigma) - \Phi(\alpha_{1,i}/\sigma)$:
/// $$\text{LOT}_{i,\text{bps}} = 2 \cdot 10^4\,\sigma_{\text{latent}}\sqrt{\frac{\pi}{2}}\left(z_{\text{zero}} + \frac{\pi}{12}z_{\text{zero}}^3\right) + 0.45\,s_{\text{ba,bps}},$$
/// $$s_{\text{liq,bps}}^{\text{CLW}} = \beta_{\text{LOT}}\,\text{LOT}_{i,\text{bps}}\left(1 + 0.06\,D_{\text{mat}}\right)\left(1 + 0.85\,L_{\text{lev}}\right).$$
pub fn evaluate_chen_lesmond_wei_lot_zero_return_bond_illiquidity<T: Real>(
    spec: &ChenLesmondWeiLotBondIlliquiditySpec<T>,
) -> ChenLesmondWeiLotBondIlliquidityAnalytics<T> {
    let one = T::one();
    let bps = T::from_f64(10_000.0);
    let sqrt_half_pi = T::from_f64(1.253_314_137_315_500_3);
    let pi_over_12 = T::from_f64(0.261_799_387_799_149_4);

    let z = spec
        .zero_return_day_proportion_z
        .max_val(T::from_f64(0.01))
        .min_val(T::from_f64(0.95));
    let sig = spec
        .latent_daily_return_volatility_sigma
        .max_val(T::from_f64(1e-4));
    let s_ba = spec.quoted_bid_ask_spread_bps.max_val(T::from_f64(0.5));
    let d_mat = spec.bond_maturity_years_d.max_val(T::from_f64(0.25));
    let lev = spec
        .issuer_leverage_ratio_l
        .max_val(T::from_f64(0.01))
        .min_val(T::from_f64(0.95));
    let beta = spec
        .lot_yield_spread_loading_beta
        .max_val(T::from_f64(0.05));

    let probit_z = sqrt_half_pi * (z + pi_over_12 * z * z * z);
    let z_band_bps = T::from_f64(2.0) * bps * sig * probit_z;
    let lot_bps = z_band_bps + T::from_f64(0.45) * s_ba;
    let yield_comp_bps =
        beta * lot_bps * (one + T::from_f64(0.06) * d_mat) * (one + T::from_f64(0.85) * lev);

    ChenLesmondWeiLotBondIlliquidityAnalytics {
        lot_roundtrip_transaction_cost_bps: lot_bps,
        zero_return_threshold_band_bps: z_band_bps,
        lot_corporate_yield_spread_component_bps: yield_comp_bps,
    }
}

#[cfg(test)]
mod chen_lesmond_wei_lot_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn chen_lesmond_wei_lot_invariants_and_ad_parity() {
        let base = ChenLesmondWeiLotBondIlliquiditySpec {
            zero_return_day_proportion_z: 0.28_f64,
            latent_daily_return_volatility_sigma: 0.0045_f64,
            quoted_bid_ask_spread_bps: 36.0_f64,
            bond_maturity_years_d: 6.5_f64,
            issuer_leverage_ratio_l: 0.42_f64,
            lot_yield_spread_loading_beta: 0.85_f64,
        };
        let r0 = evaluate_chen_lesmond_wei_lot_zero_return_bond_illiquidity(&base);
        assert!(r0.zero_return_threshold_band_bps > 0.0);
        assert!(r0.lot_roundtrip_transaction_cost_bps > r0.zero_return_threshold_band_bps);
        assert!(r0.lot_corporate_yield_spread_component_bps > 30.0);

        let x0 = [
            base.zero_return_day_proportion_z,
            base.latent_daily_return_volatility_sigma,
            base.quoted_bid_ask_spread_bps,
            base.bond_maturity_years_d,
            base.issuer_leverage_ratio_l,
            base.lot_yield_spread_loading_beta,
        ];
        let mut xd = [Dual::<6>::from_f64(0.0); 6];
        for i in 0..6 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_chen_lesmond_wei_lot_zero_return_bond_illiquidity(
            &ChenLesmondWeiLotBondIlliquiditySpec {
                zero_return_day_proportion_z: xd[0],
                latent_daily_return_volatility_sigma: xd[1],
                quoted_bid_ask_spread_bps: xd[2],
                bond_maturity_years_d: xd[3],
                issuer_leverage_ratio_l: xd[4],
                lot_yield_spread_loading_beta: xd[5],
            },
        );
        let h = 1e-6;
        for i in 0..6 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_chen_lesmond_wei_lot_zero_return_bond_illiquidity(
                &ChenLesmondWeiLotBondIlliquiditySpec {
                    zero_return_day_proportion_z: xp[0],
                    latent_daily_return_volatility_sigma: xp[1],
                    quoted_bid_ask_spread_bps: xp[2],
                    bond_maturity_years_d: xp[3],
                    issuer_leverage_ratio_l: xp[4],
                    lot_yield_spread_loading_beta: xp[5],
                },
            );
            let rm = evaluate_chen_lesmond_wei_lot_zero_return_bond_illiquidity(
                &ChenLesmondWeiLotBondIlliquiditySpec {
                    zero_return_day_proportion_z: xm[0],
                    latent_daily_return_volatility_sigma: xm[1],
                    quoted_bid_ask_spread_bps: xm[2],
                    bond_maturity_years_d: xm[3],
                    issuer_leverage_ratio_l: xm[4],
                    lot_yield_spread_loading_beta: xm[5],
                },
            );
            let fd = (rp.lot_corporate_yield_spread_component_bps
                - rm.lot_corporate_yield_spread_component_bps)
                / (2.0 * h);
            assert!(rad.lot_corporate_yield_spread_component_bps.d[i].abs() > 1e-9);
            assert!((rad.lot_corporate_yield_spread_component_bps.d[i] - fd).abs() < 1e-4);
        }
    }
}
