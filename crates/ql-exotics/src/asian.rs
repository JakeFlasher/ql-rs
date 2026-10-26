#![allow(clippy::needless_range_loop)]
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
    const X_GL8: [f64; 8] = [
        -0.960289856497536,
        -0.796666477413627,
        -0.525532409916329,
        -0.183434642495650,
        0.183434642495650,
        0.525532409916329,
        0.796666477413627,
        0.960289856497536,
    ];
    const W_GL8: [f64; 8] = [
        0.101228536290376,
        0.222381034453374,
        0.313706645877887,
        0.362683783378362,
        0.362683783378362,
        0.313706645877887,
        0.222381034453374,
        0.101228536290376,
    ];
    let n = s.fixing_times.len().max(1);
    let nf = T::from_f64(n as f64);
    let t_end = *s.fixing_times.last().unwrap_or(&T::one());
    let df = (-s.risk_free_rate * t_end).exp();
    let sig = s.vol.abs().max_val(T::from_f64(1e-8));
    let s0 = s.spot.max_val(T::from_f64(1e-8));
    let k_strike = s.strike.max_val(T::from_f64(1e-8));
    let mut mu_i = Vec::with_capacity(n);
    let mut c_i = Vec::with_capacity(n);
    let mut mu_g = T::zero();
    let mut sig_g_sq = T::zero();
    for i in 0..n {
        let ti = s
            .fixing_times
            .get(i)
            .copied()
            .unwrap_or(T::one())
            .max_val(T::from_f64(1e-8));
        let mui =
            s0.ln() + (s.risk_free_rate - s.dividend_yield - T::from_f64(0.5) * sig * sig) * ti;
        let mut ci = T::zero();
        for j in 0..n {
            let tj = s
                .fixing_times
                .get(j)
                .copied()
                .unwrap_or(T::one())
                .max_val(T::from_f64(1e-8));
            ci += sig * sig * ti.min_val(tj);
        }
        ci /= nf;
        mu_i.push(mui);
        c_i.push(ci);
        mu_g += mui / nf;
        sig_g_sq += ci / nf;
    }
    let sig_g = sig_g_sq.max_val(T::from_f64(1e-16)).sqrt();
    let mut beta_i = Vec::with_capacity(n);
    let mut v_cond_sq = Vec::with_capacity(n);
    for i in 0..n {
        let bi = c_i[i] / sig_g_sq.max_val(T::from_f64(1e-16));
        let ti = s
            .fixing_times
            .get(i)
            .copied()
            .unwrap_or(T::one())
            .max_val(T::from_f64(1e-8));
        beta_i.push(bi);
        v_cond_sq.push((sig * sig * ti - bi * bi * sig_g_sq).max_val(T::zero()));
    }
    let fg = (mu_g + T::from_f64(0.5) * sig_g_sq).exp();
    let d1_g = ((fg / k_strike).ln() + T::from_f64(0.5) * sig_g_sq) / sig_g;
    let d2_g = d1_g - sig_g;
    let c_geo = df * (fg * d1_g.norm_cdf() - k_strike * d2_g.norm_cdf()).max_val(T::zero());
    let mut lo = -12.0_f64;
    let mut hi = 12.0_f64;
    for _ in 0..32 {
        let mid = 0.5 * (lo + hi);
        let mut m_val = 0.0;
        for i in 0..n {
            m_val +=
                (mu_i[i].val() + beta_i[i].val() * sig_g.val() * mid + 0.5 * v_cond_sq[i].val())
                    .exp()
                    / n as f64;
        }
        if m_val < k_strike.val() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mut z_star = T::from_f64(0.5 * (lo + hi));
    for _ in 0..3 {
        let mut m_val = -k_strike;
        let mut dm_val = T::zero();
        for i in 0..n {
            let term =
                (mu_i[i] + beta_i[i] * sig_g * z_star + T::from_f64(0.5) * v_cond_sq[i]).exp() / nf;
            m_val += term;
            dm_val += beta_i[i] * sig_g * term;
        }
        z_star -= m_val / dm_val.max_val(T::from_f64(1e-12));
    }
    let mut sum_f_phi = T::zero();
    let mut sum_f_dspot = T::zero();
    let mut sum_f_dvega = T::zero();
    for i in 0..n {
        let ti = s
            .fixing_times
            .get(i)
            .copied()
            .unwrap_or(T::one())
            .max_val(T::from_f64(1e-8));
        let fi = s0 * ((s.risk_free_rate - s.dividend_yield) * ti).exp() / nf;
        let arg = -z_star + beta_i[i] * sig_g;
        sum_f_phi += fi * arg.norm_cdf();
        sum_f_dspot += (fi / s0) * arg.norm_cdf();
        sum_f_dvega += fi * arg.norm_pdf() * beta_i[i] * (sig_g / sig);
    }
    let i1 = df * (sum_f_phi - k_strike * (-z_star).norm_cdf()).max_val(c_geo);
    let mut i2 = T::zero();
    let mut delta_rs = T::zero();
    if n > 1 {
        let mid_z = z_star - T::from_f64(3.0);
        for q in 0..8 {
            let zq = mid_z + T::from_f64(3.0) * T::from_f64(X_GL8[q]);
            let wq = T::from_f64(3.0) * T::from_f64(W_GL8[q]) * zq.norm_pdf();
            let mut mz = T::zero();
            let mut mi_vec = Vec::with_capacity(n);
            for i in 0..n {
                let mi =
                    (mu_i[i] + beta_i[i] * sig_g * zq + T::from_f64(0.5) * v_cond_sq[i]).exp() / nf;
                mz += mi;
                mi_vec.push(mi);
            }
            let mut var_cond = T::zero();
            for i in 0..n {
                for j in 0..n {
                    let ti = s
                        .fixing_times
                        .get(i)
                        .copied()
                        .unwrap_or(T::one())
                        .max_val(T::from_f64(1e-8));
                    let tj = s
                        .fixing_times
                        .get(j)
                        .copied()
                        .unwrap_or(T::one())
                        .max_val(T::from_f64(1e-8));
                    let cov_ij = (sig * sig * ti.min_val(tj) - beta_i[i] * beta_i[j] * sig_g_sq)
                        .max_val(T::zero());
                    var_cond += mi_vec[i] * mi_vec[j] * (cov_ij.exp() - T::one());
                }
            }
            let sd_cond = var_cond.max_val(T::from_f64(1e-16)).sqrt();
            let d_tail = (mz - k_strike) / sd_cond;
            let call_cond = (mz - k_strike) * d_tail.norm_cdf() + sd_cond * d_tail.norm_pdf();
            i2 += df * wq * call_cond.max_val(T::zero());
            delta_rs += T::from_f64(0.5) * df * wq * sd_cond;
        }
    }
    let full_price = i1 + i2;
    let upper_bound = full_price + delta_rs;
    CurranRogersShiArithmeticAsianReport {
        geometric_asian_call_price: c_geo,
        curran_jensen_lower_bound_i1: i1,
        curran_full_conditioned_price: full_price,
        rogers_shi_upper_bound: upper_bound,
        rogers_shi_bound_spread: delta_rs,
        critical_conditioning_z_star: z_star,
        delta: df * sum_f_dspot,
        vega: df * sum_f_dvega,
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
mod domain_suite_39_tests {
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
        assert!(d.d[0].abs() > 1e-6 && (d.d[0] - fd).abs() < 1e-5);
    }

    #[test]
    fn suite39_invariants() {
        let r = price_curran_rogers_shi_arithmetic_asian(&make(100.0));
        assert!(r.geometric_asian_call_price <= r.curran_jensen_lower_bound_i1);
        assert!(r.curran_jensen_lower_bound_i1 <= r.curran_full_conditioned_price);
        assert!(r.curran_full_conditioned_price <= r.rogers_shi_upper_bound);
    }
}
