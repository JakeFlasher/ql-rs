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

#[inline]
fn is_breached<T: Real>(s: &ParisianBarrierOptionSpec<T>, spot: &T) -> bool {
    match s.barrier_type {
        BarrierType::UpAndOut | BarrierType::UpAndIn => spot.val() >= s.barrier.val(),
        BarrierType::DownAndOut | BarrierType::DownAndIn => spot.val() <= s.barrier.val(),
    }
}

#[inline]
fn is_knock_in<T: Real>(s: &ParisianBarrierOptionSpec<T>) -> bool {
    matches!(
        s.barrier_type,
        BarrierType::UpAndIn | BarrierType::DownAndIn
    )
}

#[inline]
fn vanilla_payoff<T: Real>(s: &ParisianBarrierOptionSpec<T>, spot: T) -> T {
    if s.option_type.is_call() {
        (spot - s.strike).max_val(T::zero())
    } else {
        (s.strike - spot).max_val(T::zero())
    }
}

/// Recombining binomial lattice with an additional discrete excursion-clock
/// state.  The clock is measured in observation steps; consecutive clocks
/// reset on an out-of-barrier observation, while cumulative clocks retain all
/// time spent outside.  A zero window therefore reduces to an ordinary
/// discretely-monitored barrier option.
pub fn price_parisian_option_lattice<T: Real>(
    s: &ParisianBarrierOptionSpec<T>,
    n_steps: usize,
) -> T {
    let n = n_steps.max(16);
    if !s.expiry.val().is_finite() || s.expiry.val() < 0.0 {
        return T::zero();
    }
    if s.expiry.val() <= 0.0 {
        let pay = vanilla_payoff(s, s.spot);
        if s.window_duration.val() <= 0.0 && is_breached(s, &s.spot) {
            return if is_knock_in(s) { pay } else { s.rebate };
        }
        return if is_knock_in(s) { T::zero() } else { pay };
    }
    // A zero-volatility lattice has u == d and an undefined risk-neutral
    // probability.  Use the deterministic terminal payoff with the barrier
    // state sampled at each step instead of manufacturing a NaN.
    if s.vol.val().abs() < 1e-14 {
        let dt = s.expiry / n as f64;
        let mut x = s.spot;
        let mut clock = 0usize;
        let zero_window = s.window_duration.val() <= 0.0;
        // A zero excursion window means a breach at an observation is
        // sufficient.  Represent it with a one-step clock (rather than a
        // zero threshold, for which every state would be absorbing).
        let window_steps = if zero_window {
            1
        } else {
            ((s.window_duration.val() / dt.val()).ceil() as usize).min(n)
        };
        if zero_window && is_breached(s, &x) {
            return if is_knock_in(s) {
                vanilla_payoff(s, s.spot)
            } else {
                s.rebate
            };
        }
        let mut knocked = false;
        for _ in 0..n {
            x *= ((s.rate - s.div_yield) * dt).exp();
            let b = is_breached(s, &x);
            clock = match s.clock_kind {
                ExcursionClockKind::ParisianConsecutive => {
                    if b {
                        clock.saturating_add(1)
                    } else {
                        0
                    }
                }
                ExcursionClockKind::ParasianCumulative => {
                    if b {
                        clock.saturating_add(1)
                    } else {
                        clock
                    }
                }
            };
            knocked |= clock >= window_steps;
        }
        let payoff = vanilla_payoff(s, x);
        let terminal = if is_knock_in(s) {
            if knocked { payoff } else { T::zero() }
        } else if knocked {
            s.rebate
        } else {
            payoff
        };
        return (-s.rate * s.expiry).exp() * terminal;
    }
    let dt = s.expiry / T::from_f64(n as f64);
    let u = (s.vol * dt.sqrt()).exp();
    let d = T::one() / u;
    let a = ((s.rate - s.div_yield) * dt).exp();
    let p = ((a - d) / (u - d)).max_val(T::zero()).min_val(T::one());
    let zero_window = s.window_duration.val() <= 0.0;
    let window_steps = if zero_window {
        // See the zero-volatility branch above.  The first breach transition
        // moves the clock from 0 to 1 and therefore triggers immediately.
        1
    } else {
        ((s.window_duration.val() / dt.val()).ceil() as usize).min(n)
    };
    let clocks = window_steps.saturating_add(1);
    // Values at maturity and the preceding time layer.  The flat arrays are
    // indexed [spot-node][clock-state] and avoid a cubic path enumeration.
    let mut next = vec![T::zero(); (n + 1) * clocks];
    for j in 0..=n {
        let st = s.spot * u.powf((2 * j) as f64 - n as f64);
        let pay = vanilla_payoff(s, st);
        for c in 0..clocks {
            let knocked = c >= window_steps;
            next[j * clocks + c] = if is_knock_in(s) {
                if knocked { pay } else { T::zero() }
            } else if knocked {
                s.rebate
            } else {
                pay
            };
        }
    }
    for step in (0..n).rev() {
        let mut cur = vec![T::zero(); (step + 1) * clocks];
        for j in 0..=step {
            let st_up = s.spot * u.powf((2 * (j + 1)) as f64 - (step + 1) as f64);
            let st_dn = s.spot * u.powf((2 * j) as f64 - (step + 1) as f64);
            let bu = is_breached(s, &st_up);
            let bd = is_breached(s, &st_dn);
            for c in 0..clocks {
                let cu = match s.clock_kind {
                    ExcursionClockKind::ParisianConsecutive => {
                        // The threshold state is absorbing once the
                        // Parisian excursion has completed.
                        if c >= window_steps {
                            window_steps
                        } else if bu {
                            c.saturating_add(1).min(window_steps)
                        } else {
                            0
                        }
                    }
                    ExcursionClockKind::ParasianCumulative => {
                        if bu {
                            c.saturating_add(1).min(window_steps)
                        } else {
                            c
                        }
                    }
                };
                let cd = match s.clock_kind {
                    ExcursionClockKind::ParisianConsecutive => {
                        if c >= window_steps {
                            window_steps
                        } else if bd {
                            c.saturating_add(1).min(window_steps)
                        } else {
                            0
                        }
                    }
                    ExcursionClockKind::ParasianCumulative => {
                        if bd {
                            c.saturating_add(1).min(window_steps)
                        } else {
                            c
                        }
                    }
                };
                cur[j * clocks + c] = (-s.rate * dt).exp()
                    * (p * next[(j + 1) * clocks + cu] + (T::one() - p) * next[j * clocks + cd]);
            }
        }
        next = cur;
    }
    // A zero-duration contract is nevertheless immediately knocked in/out if
    // the spot starts on the breached side.  For a positive window the clock
    // starts at zero and accrues only subsequent observations.
    if zero_window && is_breached(s, &s.spot) {
        if is_knock_in(s) {
            return vanilla_payoff(s, s.spot);
        }
        return s.rebate;
    }
    next[0]
}

/// Antithetic Monte Carlo estimator with a smooth barrier indicator.  The
/// random variates are Box--Muller normals (the old implementation used a
/// uniform draw, materially understating diffusion variance).
pub fn price_parisian_option_mc_ad<T: Real>(
    s: &ParisianBarrierOptionSpec<T>,
    n_paths: usize,
    n_steps: usize,
    smoothing_eps: f64,
    seed: u64,
) -> T {
    if !s.expiry.val().is_finite() || s.expiry.val() < 0.0 {
        return T::zero();
    }
    if s.expiry.val() <= 0.0 {
        let pay = vanilla_payoff(s, s.spot);
        if s.window_duration.val() <= 0.0 && is_breached(s, &s.spot) {
            return if is_knock_in(s) { pay } else { s.rebate };
        }
        return if is_knock_in(s) { T::zero() } else { pay };
    }
    let n = n_steps.max(1);
    let paths = n_paths.max(2);
    let dt = s.expiry / T::from_f64(n as f64);
    let mut total = T::zero();
    for pth in 0..paths {
        let sign = if pth % 2 == 0 { 1.0 } else { -1.0 };
        let mut seedv = seed.wrapping_add((pth / 2) as u64);
        let mut x = s.spot;
        let mut clock = T::zero();
        let mut knocked = if s.window_duration.val() <= 0.0 {
            T::from_f64(if is_breached(s, &x) { 1.0 } else { 0.0 })
        } else {
            T::zero()
        };
        for _ in 0..n {
            seedv = seedv.wrapping_mul(6364136223846793005).wrapping_add(1);
            let u1 = (((seedv >> 11) as f64) / 9007199254740992.0).max(1e-15);
            seedv = seedv.wrapping_mul(6364136223846793005).wrapping_add(1);
            let u2 = (seedv >> 11) as f64 / 9007199254740992.0;
            let z = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos() * sign;
            x *= ((s.rate - s.div_yield - s.vol * s.vol * 0.5) * dt + s.vol * dt.sqrt() * z).exp();
            let breach = match s.barrier_type {
                BarrierType::UpAndOut | BarrierType::UpAndIn => {
                    (x - s.barrier).smooth_heaviside_c2(T::zero(), smoothing_eps.max(1e-8))
                }
                BarrierType::DownAndOut | BarrierType::DownAndIn => {
                    (s.barrier - x).smooth_heaviside_c2(T::zero(), smoothing_eps.max(1e-8))
                }
            };
            clock = match s.clock_kind {
                ExcursionClockKind::ParisianConsecutive => breach * (clock + dt),
                ExcursionClockKind::ParasianCumulative => clock + breach * dt,
            };
            if s.window_duration.val() <= 0.0 {
                // For an instantaneous excursion window, a breach at the
                // observation itself is the stopping event.  Applying the
                // smooth clock gate at threshold zero would otherwise assign
                // an erroneous 0.5 weight to paths that never breach.
                knocked = knocked.max_val(breach);
            } else {
                let clock_hw = (s.expiry.val() / n as f64).max(smoothing_eps).max(1e-8);
                let gate = clock
                    .smooth_heaviside_c2(s.window_duration, clock_hw)
                    .max_val(T::zero());
                knocked = knocked.max_val(gate);
            }
        }
        let pay = vanilla_payoff(s, x);
        let terminal = if is_knock_in(s) {
            knocked * pay
        } else {
            knocked * s.rebate + (T::one() - knocked) * pay
        };
        total += (-s.rate * s.expiry).exp() * terminal;
    }
    total / paths as f64
}
