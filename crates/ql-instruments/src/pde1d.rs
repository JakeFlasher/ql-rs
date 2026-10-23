use crate::options::{AmericanOption, BarrierType, EuropeanOption, OptionType};
use ql_ad::Real;
use ql_math::solve_tridiagonal;
#[derive(Debug, Clone, Copy)]
pub struct Fd1dBlackScholesEngine {
    pub grid_points: usize,
    pub time_steps: usize,
}
impl Default for Fd1dBlackScholesEngine {
    fn default() -> Self {
        Self {
            grid_points: 200,
            time_steps: 200,
        }
    }
}
impl Fd1dBlackScholesEngine {
    pub fn new(grid_points: usize, time_steps: usize) -> Self {
        Self {
            grid_points: grid_points.max(40),
            time_steps: time_steps.max(2),
        }
    }
    #[allow(clippy::too_many_arguments)]
    pub fn price<T: Real>(
        &self,
        spot: T,
        strike: T,
        expiry: T,
        vol: T,
        rate: T,
        dividend: T,
        kind: OptionType,
    ) -> T {
        self.fd(spot, strike, expiry, vol, rate, dividend, kind, None, false)
    }
    #[allow(
        clippy::too_many_arguments,
        clippy::needless_range_loop,
        clippy::collapsible_if
    )]
    fn fd<T: Real>(
        &self,
        spot: T,
        strike: T,
        expiry: T,
        vol: T,
        rate: T,
        dividend: T,
        kind: OptionType,
        barrier: Option<(T, BarrierType)>,
        american: bool,
    ) -> T {
        let n = self.grid_points.max(40);
        let m = self.time_steps.max(2);
        let s0 = spot.val().max(1e-8);
        let k = strike.val().max(1e-8);
        let smin = (k * 0.01).min(s0 * 0.1).max(1e-6);
        let smax = (4.0 * s0).max(4.0 * k).max(
            barrier
                .map(|(b, t)| {
                    if matches!(t, BarrierType::UpAndOut | BarrierType::UpAndIn) {
                        2.0 * b.val()
                    } else {
                        4.0 * b.val()
                    }
                })
                .unwrap_or(0.0),
        );
        let mut xmin = smin.ln();
        let dx = (smax.ln() - xmin) / n as f64;
        // Translate the log grid so a monitored barrier is represented by an
        // exact node rather than an O(dx) quantized level.
        if let Some((barrier_level, _)) = barrier {
            let j = ((barrier_level.val().ln() - xmin) / dx).round();
            xmin = barrier_level.val().ln() - j * dx;
        }
        let mut v = vec![T::zero(); n + 1];
        for j in 0..=n {
            let sj = T::from_f64((xmin + j as f64 * dx).exp());
            v[j] = if kind.is_call() {
                (sj - strike).max_val(T::zero())
            } else {
                (strike - sj).max_val(T::zero())
            };
            if let Some((b, bt)) = barrier
                && matches!(bt, BarrierType::UpAndOut | BarrierType::DownAndOut)
                && ((matches!(bt, BarrierType::UpAndOut) && sj.val() >= b.val())
                    || (matches!(bt, BarrierType::DownAndOut) && sj.val() <= b.val()))
            {
                v[j] = T::zero();
            }
        }
        let mut lower = vec![T::zero(); n - 1];
        let mut diag = vec![T::zero(); n - 1];
        let mut upper = vec![T::zero(); n - 1];
        let mut rhs = vec![T::zero(); n - 1];
        let mut out = vec![T::zero(); n - 1];
        let mut scratch = vec![T::zero(); n - 1];
        let mut old = vec![T::zero(); n + 1];
        let dt = expiry / T::from_f64(m as f64);
        let mut tau = T::zero();
        for step in 0..=m {
            let (h, theta) = if step < 2 {
                (dt * 0.5, T::one())
            } else {
                (dt, T::from_f64(0.5))
            };
            old.copy_from_slice(&v);
            let sig2 = vol * vol;
            let drift = rate - dividend - T::from_f64(0.5) * sig2;
            for j in 1..n {
                let aa =
                    T::from_f64(0.5) * sig2 / T::from_f64(dx * dx) - drift / T::from_f64(2.0 * dx);
                let cc =
                    T::from_f64(0.5) * sig2 / T::from_f64(dx * dx) + drift / T::from_f64(2.0 * dx);
                let bb = -(aa + cc) - rate;
                let th = theta * h;
                lower[j - 1] = -th * aa;
                diag[j - 1] = T::one() - th * bb;
                upper[j - 1] = -th * cc;
                let ex = (T::one() - theta) * h;
                rhs[j - 1] = old[j] + ex * (aa * old[j - 1] + bb * old[j] + cc * old[j + 1]);
            }
            let tnew = (tau + h).min_val(expiry);
            let callb =
                T::from_f64(smax) * (-dividend * tnew).exp() - strike * (-rate * tnew).exp();
            // At the lower spot boundary an American put is exercised
            // immediately, hence its exact asymptotic value is K-S (not the
            // discounted European K e^{-r tau}-S e^{-q tau}).  Using the
            // European boundary here suppresses early exercise near S=0 and
            // biases the American price/delta on coarse grids.
            let putb = if american {
                (strike - T::from_f64(smin)).max_val(T::zero())
            } else {
                strike * (-rate * tnew).exp()
            };
            rhs[0] -= lower[0] * if kind.is_call() { T::zero() } else { putb };
            rhs[n - 2] -= upper[n - 2] * if kind.is_call() { callb } else { T::zero() };
            if let Some((barrier, bt)) = barrier {
                for j in 1..n {
                    let sj = (xmin + j as f64 * dx).exp();
                    let knocked = (matches!(bt, BarrierType::UpAndOut) && sj >= barrier.val())
                        || (matches!(bt, BarrierType::DownAndOut) && sj <= barrier.val());
                    if knocked {
                        lower[j - 1] = T::zero();
                        diag[j - 1] = T::one();
                        upper[j - 1] = T::zero();
                        rhs[j - 1] = T::zero();
                    }
                }
            }
            solve_tridiagonal(&lower, &diag, &upper, &rhs, &mut out, &mut scratch);
            v[1..n].copy_from_slice(&out);
            v[0] = if kind.is_call() { T::zero() } else { putb };
            v[n] = if kind.is_call() { callb } else { T::zero() };
            if american {
                for j in 1..n {
                    let sj = T::from_f64((xmin + j as f64 * dx).exp());
                    let intr = if kind.is_call() {
                        (sj - strike).max_val(T::zero())
                    } else {
                        (strike - sj).max_val(T::zero())
                    };
                    v[j] = v[j].max_val(intr);
                }
            }
            if let Some((b, bt)) = barrier
                && matches!(bt, BarrierType::UpAndOut | BarrierType::DownAndOut)
            {
                for j in 0..=n {
                    let sj = (xmin + j as f64 * dx).exp();
                    if (matches!(bt, BarrierType::UpAndOut) && sj >= b.val())
                        || (matches!(bt, BarrierType::DownAndOut) && sj <= b.val())
                    {
                        v[j] = T::zero();
                    }
                }
            }
            tau = (tau + h).min_val(expiry);
        }
        let x = spot.ln();
        let z = ((x.val() - xmin) / dx).clamp(1.0, n as f64 - 1.0);
        let i = z.floor() as usize;
        let xm = T::from_f64(xmin + (i - 1) as f64 * dx);
        let x0 = T::from_f64(xmin + i as f64 * dx);
        let xp = T::from_f64(xmin + (i + 1) as f64 * dx);
        let lm = (x - x0) * (x - xp) / T::from_f64(2.0 * dx * dx);
        let l0 = -(x - xm) * (x - xp) / T::from_f64(dx * dx);
        let lp = (x - xm) * (x - x0) / T::from_f64(2.0 * dx * dx);
        lm * v[i - 1] + l0 * v[i] + lp * v[i + 1]
    }
    #[allow(clippy::too_many_arguments)]
    pub fn price_barrier<T: Real>(
        &self,
        spot: T,
        strike: T,
        barrier: T,
        expiry: T,
        vol: T,
        rate: T,
        kind: OptionType,
        bt: BarrierType,
    ) -> T {
        if matches!(bt, BarrierType::UpAndIn | BarrierType::DownAndIn) {
            let vanilla = self.price(spot, strike, expiry, vol, rate, T::zero(), kind);
            let out = self.fd(
                spot,
                strike,
                expiry,
                vol,
                rate,
                T::zero(),
                kind,
                Some((
                    barrier,
                    if matches!(bt, BarrierType::UpAndIn) {
                        BarrierType::UpAndOut
                    } else {
                        BarrierType::DownAndOut
                    },
                )),
                false,
            );
            vanilla - out
        } else {
            self.fd(
                spot,
                strike,
                expiry,
                vol,
                rate,
                T::zero(),
                kind,
                Some((barrier, bt)),
                false,
            )
        }
    }
    pub fn price_european<T: Real>(&self, o: &EuropeanOption) -> T {
        // `EuropeanOption` is parameterized by a *forward* and a discount
        // factor (the same convention as Black--76), whereas `price` solves
        // the spot Black--Scholes PDE.  Run the PDE in the zero-rate forward
        // measure and apply the supplied discount factor once.  Passing the
        // forward straight through without this factor returned an
        // undiscounted value and silently disagreed with `EuropeanOption::price`.
        T::from_f64(o.discount)
            * self.price(
                T::from_f64(o.forward),
                T::from_f64(o.strike),
                T::from_f64(o.expiry),
                T::from_f64(o.volatility),
                T::zero(),
                T::zero(),
                o.option_type,
            )
    }
    pub fn price_american<T: Real>(&self, o: &AmericanOption) -> T {
        self.fd(
            T::from_f64(o.spot),
            T::from_f64(o.strike),
            T::from_f64(o.expiry),
            T::from_f64(o.volatility),
            T::from_f64(o.rate),
            T::zero(),
            o.option_type,
            None,
            true,
        )
    }
    pub fn price_american_ad<T: Real>(
        &self,
        spot: T,
        strike: T,
        expiry: T,
        vol: T,
        rate: T,
        kind: OptionType,
    ) -> T {
        self.fd(spot, strike, expiry, vol, rate, T::zero(), kind, None, true)
    }
    pub fn price_barrier_option<T: Real>(&self, o: &crate::options::BarrierOption) -> T {
        self.price_barrier(
            T::from_f64(o.underlying),
            T::from_f64(o.strike),
            T::from_f64(o.barrier),
            T::from_f64(o.expiry),
            T::from_f64(o.volatility),
            T::from_f64(o.rate),
            o.option_type,
            o.barrier_type,
        )
    }
}
