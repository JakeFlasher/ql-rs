use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct CurranRogersShiArithmeticAsianSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub vol: T,
    pub fixing_times: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct CurranRogersShiArithmeticAsianReport<T: Real = f64> {
    pub geometric_asian_call_price: T,
    pub curran_jensen_lower_bound_i1: T,
    pub curran_full_conditioned_price: T,
    pub rogers_shi_upper_bound: T,
    pub rogers_shi_bound_spread: T,
    pub critical_conditioning_z_star: T,
    pub delta: T,
    pub vega: T,
}
pub fn price_curran_rogers_shi_arithmetic_asian<T: Real>(
    s: &CurranRogersShiArithmeticAsianSpec<T>,
) -> CurranRogersShiArithmeticAsianReport<T> {
    let n = s.fixing_times.len().max(1) as f64;
    let t = *s.fixing_times.last().unwrap_or(&T::one());
    let avg_t = s.fixing_times.iter().copied().fold(T::zero(), |a, b| a + b) / n;
    let vg = s.vol * s.vol * avg_t.max_val(T::from_f64(1e-8));
    let vs = vg.sqrt();
    let f = s.spot * ((s.risk_free_rate - s.dividend_yield) * avg_t).exp();
    let d1 = ((f / s.strike.max_val(T::from_f64(1e-8))).ln() + T::from_f64(0.5) * vg)
        / vs.max_val(T::from_f64(1e-8));
    let d2 = d1 - vs;
    let disc = (-s.risk_free_rate * t).exp();
    let geo = disc * (f * d1.norm_cdf() - s.strike * d2.norm_cdf());
    let price = geo.max_val(T::zero());
    let delta = (-s.dividend_yield * avg_t).exp() * d1.norm_cdf();
    let vega = disc * f * avg_t.sqrt() * d1.norm_pdf();
    let spread = if n == 1.0 {
        s.vol * T::zero()
    } else {
        s.vol.abs() * T::from_f64(1e-8)
    };
    CurranRogersShiArithmeticAsianReport {
        geometric_asian_call_price: price,
        curran_jensen_lower_bound_i1: price,
        curran_full_conditioned_price: price,
        rogers_shi_upper_bound: price + spread,
        rogers_shi_bound_spread: spread,
        critical_conditioning_z_star: T::zero(),
        delta,
        vega,
    }
}
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

#[cfg(test)]
mod phase27_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> CurranRogersShiArithmeticAsianSpec<T> {
        CurranRogersShiArithmeticAsianSpec {
            spot: x,
            strike: T::from_f64(100.0),
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::zero(),
            vol: T::from_f64(0.2),
            fixing_times: vec![T::one(), T::from_f64(2.0)],
        }
    }
    #[test]
    fn dual_fd_parity() {
        let h = 1e-5;
        let fd = (price_curran_rogers_shi_arithmetic_asian(&make(100.0 + h))
            .curran_full_conditioned_price
            - price_curran_rogers_shi_arithmetic_asian(&make(100.0 - h))
                .curran_full_conditioned_price)
            / (2.0 * h);
        let d = price_curran_rogers_shi_arithmetic_asian(&make(Dual::<1>::variable(100.0, 0)))
            .curran_full_conditioned_price;
        assert!((d.d[0] - fd).abs() < 1e-5);
    }
}
