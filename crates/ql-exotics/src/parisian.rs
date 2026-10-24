use ql_ad::Real;
use ql_instruments::{BarrierType, OptionType};
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ExcursionClockKind {
    ParisianConsecutive,
    ParasianCumulative,
}
#[derive(Copy, Clone, Debug)]
pub struct ParisianBarrierOptionSpec<T: Real> {
    pub option_type: OptionType,
    pub barrier_type: BarrierType,
    pub clock_kind: ExcursionClockKind,
    pub spot: T,
    pub strike: T,
    pub barrier: T,
    pub rebate: T,
    pub window_duration: T,
    pub rate: T,
    pub div_yield: T,
    pub vol: T,
    pub expiry: T,
}
pub fn price_parisian_option_lattice<T: Real>(
    s: &ParisianBarrierOptionSpec<T>,
    n_steps: usize,
) -> T {
    let n = n_steps.max(16);
    let dt = s.expiry / T::from_f64(n as f64);
    let u = (s.vol * dt.sqrt()).exp();
    let d = T::one() / u;
    let a = ((s.rate - s.div_yield) * dt).exp();
    let p = (a - d) / (u - d);
    let mut v = vec![T::zero(); n + 1];
    for (j, cell) in v.iter_mut().enumerate() {
        let st = s.spot * u.powf((2 * j) as f64 - n as f64);
        *cell = if s.option_type.is_call() {
            (st - s.strike).max_val(T::zero())
        } else {
            (s.strike - st).max_val(T::zero())
        };
    }
    for step in (0..n).rev() {
        for j in 0..=step {
            let st = s.spot * u.powf((2 * j) as f64 - step as f64);
            let breach = match s.barrier_type {
                BarrierType::UpAndOut | BarrierType::UpAndIn => st.val() >= s.barrier.val(),
                BarrierType::DownAndOut | BarrierType::DownAndIn => st.val() <= s.barrier.val(),
            };
            let mut x = (p * v[j + 1] + (T::one() - p) * v[j]) * (-s.rate * dt).exp();
            if breach
                && matches!(
                    s.barrier_type,
                    BarrierType::UpAndOut | BarrierType::DownAndOut
                )
            {
                x = s.rebate;
            }
            v[j] = x;
        }
    }
    v[0]
}
pub fn price_parisian_option_mc_ad<T: Real>(
    s: &ParisianBarrierOptionSpec<T>,
    n_paths: usize,
    n_steps: usize,
    smoothing_eps: f64,
    seed: u64,
) -> T {
    let n = n_steps.max(1);
    let dt = s.expiry / T::from_f64(n as f64);
    let mut total = T::zero();
    for pth in 0..n_paths.max(2) {
        let sign = if pth % 2 == 0 { 1.0 } else { -1.0 };
        let mut seedv = seed.wrapping_add((pth / 2) as u64);
        let mut x = s.spot;
        let mut clock = T::zero();
        for _ in 0..n {
            seedv = seedv.wrapping_mul(6364136223846793005).wrapping_add(1);
            let z = (((seedv >> 11) as f64) / 9007199254740992.0 - 0.5) * 2.0 * sign;
            x *= ((s.rate - s.div_yield - s.vol * s.vol * 0.5) * dt + s.vol * dt.sqrt() * z).exp();
            let w = T::from_f64(smoothing_eps.max(1e-3)) * s.barrier.abs();
            let breach = match s.barrier_type {
                BarrierType::UpAndOut | BarrierType::UpAndIn => {
                    (x - s.barrier + w).smooth_heaviside_c2(T::zero(), 1.0)
                }
                BarrierType::DownAndOut | BarrierType::DownAndIn => {
                    (s.barrier - x + w).smooth_heaviside_c2(T::zero(), 1.0)
                }
            };
            clock = match s.clock_kind {
                ExcursionClockKind::ParisianConsecutive => breach * (clock + dt),
                ExcursionClockKind::ParasianCumulative => clock + breach * dt,
            };
        }
        let gate =
            clock.smooth_heaviside_c2(s.window_duration, s.window_duration.val().abs().max(1e-6));
        let knocked = match s.barrier_type {
            BarrierType::UpAndOut | BarrierType::DownAndOut => gate,
            BarrierType::UpAndIn | BarrierType::DownAndIn => T::one() - gate,
        };
        let pay = if s.option_type.is_call() {
            (x - s.strike).max_val(T::zero())
        } else {
            (s.strike - x).max_val(T::zero())
        };
        total += (-s.rate * s.expiry).exp() * (knocked * pay + (T::one() - knocked) * s.rebate);
    }
    total / n_paths.max(2) as f64
}
