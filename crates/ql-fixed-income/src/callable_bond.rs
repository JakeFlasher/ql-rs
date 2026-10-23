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
        self.schedule.iter().fold(p0, |v, c| {
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
            if c.is_call { v - option } else { v + option }
        })
    }
    pub fn tree_callable_bond_price(&self, short_rate: T, rate_vol: T, oas: T, steps: usize) -> T {
        let n = steps.max(1);
        let dt = self.maturity / n as f64;
        let vol_step = rate_vol * dt.sqrt();
        let mut values = vec![T::zero(); n + 1];
        // Principal is paid at maturity; the final coupon is added at the
        // maturity node just like every other coupon below.
        values.fill(self.face_amount);
        for (idx, &t) in self.coupon_times.iter().enumerate() {
            if (t - self.maturity).val().abs() < dt.val().abs() * 0.5 + 1e-12 {
                let cf = self.cashflow_at(idx) - self.face_amount;
                values.fill(self.face_amount + cf);
            }
        }
        let exercise_at = |j: usize| -> Vec<&BondCallability<T>> {
            let time = dt * j as f64;
            self.schedule
                .iter()
                .filter(|c| (c.exercise_time - time).val().abs() <= dt.val().abs() * 0.5 + 1e-12)
                .collect()
        };
        for j in (0..n).rev() {
            let time = dt * j as f64;
            let mut next = vec![T::zero(); j + 1];
            for i in 0..=j {
                let short_up = short_rate + oas + vol_step * (2.0 * (i + 1) as f64 - j as f64);
                let short_dn = short_rate + oas + vol_step * (2.0 * i as f64 - j as f64);
                let disc_up = (-short_up * dt).exp();
                let disc_dn = (-short_dn * dt).exp();
                next[i] = (values[i] * disc_dn + values[i + 1] * disc_up) * 0.5;
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
                if (ct - time).val().abs() <= dt.val().abs() * 0.5 + 1e-12 {
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
