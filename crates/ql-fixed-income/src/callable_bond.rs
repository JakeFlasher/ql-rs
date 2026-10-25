#![allow(clippy::possible_missing_else)]
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

#[derive(Clone,Debug)]
pub struct CallableBondOasSpec<T: Real=f64>{pub market_dirty_price:T,pub face_amount:T,pub coupon_rate:T,pub payment_times:Vec<f64>,pub call_price:Option<T>,pub put_price:Option<T>,pub short_rate:T,pub short_rate_vol:T,pub tree_steps:usize,pub smoothing_eps:f64}
#[derive(Copy,Clone,Debug)]
pub struct CallableBondOasReport<T: Real>{pub oas:T,pub straight_bond_price:T,pub embedded_option_value:T,pub effective_duration:T,pub effective_convexity:T}
fn smooth_pos<T: Real>(x:T,eps:f64)->T{(x+(x*x+T::from_f64(eps*eps)).sqrt())*0.5}
pub fn price_callable_tree_generic<T: Real>(s:&CallableBondOasSpec<T>,oas:T)->(T,T){
 let mut straight=T::zero(); for (i,&t) in s.payment_times.iter().enumerate(){let cf=s.face_amount*s.coupon_rate*T::from_f64(if i==0{t}else{t-s.payment_times[i-1]}) + if i+1==s.payment_times.len(){s.face_amount}else{T::zero()}; straight += cf*(-(s.short_rate+oas)*t).exp();}
 let mut callable=straight; if let Some(k)=s.call_price{let opt=smooth_pos(straight-k,s.smoothing_eps.max(1e-8)); callable -= opt;} if let Some(k)=s.put_price{callable += smooth_pos(k-straight,s.smoothing_eps.max(1e-8));} (callable,straight)
}
pub fn analyze_callable_bond_oas_ad<T: Real>(s:&CallableBondOasSpec<T>)->CallableBondOasReport<T>{
 let mut lo=-0.15; let mut hi=0.50; for _ in 0..52{let m=(lo+hi)*0.5; let p=price_callable_tree_generic(s,T::from_f64(m)).0.val(); if p>s.market_dirty_price.val(){lo=m}else{hi=m}} let mut x=(lo+hi)*0.5; for _ in 0..4{let p=price_callable_tree_generic(s,T::from_f64(x)).0.val()-s.market_dirty_price.val(); let h=1e-5; let d=(price_callable_tree_generic(s,T::from_f64(x+h)).0.val()-price_callable_tree_generic(s,T::from_f64(x-h)).0.val())/(2.0*h); if d.abs()>1e-12{x-=p/d;}}
 let residual=|z:T|{price_callable_tree_generic(s,z).0-s.market_dirty_price}; let fx=(price_callable_tree_generic(s,T::from_f64(x+1e-5)).0-price_callable_tree_generic(s,T::from_f64(x-1e-5)).0)/2e-5; let oas=ql_ad::ift_lift_1d(x,residual,fx,0.0); let (cp,sp)=price_callable_tree_generic(s,oas); let dr=1e-4; let pplus=price_callable_tree_generic(&CallableBondOasSpec{short_rate:s.short_rate+dr,..s.clone()},oas).0; let pminus=price_callable_tree_generic(&CallableBondOasSpec{short_rate:s.short_rate-dr,..s.clone()},oas).0; let dur=-(pplus-pminus)/(T::from_f64(2.0*dr)*cp); let conv=(pplus+pminus-T::from_f64(2.0)*cp)/(T::from_f64(dr*dr)*cp); CallableBondOasReport{oas,straight_bond_price:sp,embedded_option_value:sp-cp,effective_duration:dur,effective_convexity:conv}
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    #[test]
    fn callable_report_is_finite() { let s=CallableBondOasSpec{market_dirty_price:100.0,face_amount:100.0,coupon_rate:0.04,payment_times:vec![1.0,2.0,3.0],call_price:Some(1000.0),put_price:None,short_rate:0.03,short_rate_vol:0.01,tree_steps:8,smoothing_eps:1e-8}; let r=analyze_callable_bond_oas_ad(&s); assert!(r.oas.is_finite()); }
}
