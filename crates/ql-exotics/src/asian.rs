use ql_ad::Real;
use ql_instruments::OptionType;
use ql_vol::black76_price;

fn call_flag(option_type: OptionType) -> bool {
    option_type.is_call()
}

pub fn price_geometric_asian_continuous<T: Real>(
    spot: T,
    strike: T,
    rate: T,
    div_yield: T,
    vol: T,
    expiry: T,
    option_type: OptionType,
) -> T {
    let sigma_g = vol / 3.0_f64.sqrt();
    let b_g = (rate - div_yield - vol * vol / 6.0) * 0.5;
    let forward = spot * (b_g * expiry).exp();
    black76_price(
        forward,
        strike,
        sigma_g,
        expiry,
        (-rate * expiry).exp(),
        call_flag(option_type),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn price_geometric_asian_discrete<T: Real>(
    spot: T,
    strike: T,
    rate: T,
    div_yield: T,
    vol: T,
    fixing_times: &[T],
    expiry: T,
    option_type: OptionType,
) -> T {
    let n = fixing_times.len();
    if n == 0 {
        return T::zero();
    }
    let nf = n as f64;
    let b = rate - div_yield - vol * vol * 0.5;
    let sum_t = fixing_times.iter().copied().fold(T::zero(), |a, t| a + t);
    let mu = spot.ln() + b * (sum_t / nf);
    let mut var = T::zero();
    for (i, t) in fixing_times.iter().enumerate() {
        var += T::from_f64(2.0 * (n - i - 1) as f64 + 1.0) * *t;
    }
    var *= vol * vol / (nf * nf);
    let forward = (mu + var * 0.5).exp();
    let vol_eff = if expiry.val() > 0.0 {
        (var / expiry).sqrt()
    } else {
        T::zero()
    };
    black76_price(
        forward,
        strike,
        vol_eff,
        expiry,
        (-rate * expiry).exp(),
        call_flag(option_type),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn price_turnbull_wakeman_arithmetic_asian<T: Real>(
    spot: T,
    strike: T,
    rate: T,
    div_yield: T,
    vol: T,
    future_fixing_times: &[T],
    expiry: T,
    past_fixings_count: usize,
    running_average: T,
    option_type: OptionType,
) -> T {
    let n = future_fixing_times.len();
    if n == 0 {
        return T::zero();
    }
    let m = past_fixings_count as f64;
    let total = m + n as f64;
    let b = rate - div_yield;
    let m1 = future_fixing_times
        .iter()
        .copied()
        .map(|t| (b * t).exp())
        .fold(T::zero(), |a, x| a + x)
        * spot
        / n as f64;
    let mut m2 = T::zero();
    for (i, ti) in future_fixing_times.iter().enumerate() {
        m2 += ((b * 2.0 + vol * vol) * *ti).exp();
        for tj in future_fixing_times.iter().take(i) {
            m2 += (b * (*ti + *tj) + vol * vol * *tj).exp() * 2.0;
        }
    }
    m2 *= spot * spot / (n as f64 * n as f64);
    let variance = (m2 / (m1 * m1)).ln().max_val(T::zero());
    let vol_eff = if expiry.val() > 0.0 {
        (variance / expiry).sqrt()
    } else {
        T::zero()
    };
    let k_eff =
        T::from_f64(total / n as f64) * strike - T::from_f64(m / n as f64) * running_average;
    let scale = T::from_f64(n as f64 / total);
    let df = (-rate * expiry).exp();
    if k_eff.val() <= 0.0 {
        let intrinsic = if option_type.is_call() {
            (m1 - k_eff).max_val(T::zero())
        } else {
            T::zero()
        };
        return scale * df * intrinsic;
    }
    scale * black76_price(m1, k_eff, vol_eff, expiry, df, option_type.is_call())
}
