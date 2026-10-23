use crate::LevyModel;
use ql_ad::Real;
use ql_instruments::OptionType;
fn coeff_put<T: Real>(a: T, w: T, n: usize) -> T {
    let om = T::from_f64(n as f64 * std::f64::consts::PI) / w;
    let psi = if n == 0 { -a } else { (-om * a).sin() / om };
    let chi = ((-om * a).cos() - a.exp() + om * (-om * a).sin()) / (T::one() + om * om);
    T::from_f64(2.0) / w * (-chi + psi)
}
#[allow(clippy::too_many_arguments)]
fn price_one<T: Real, M: LevyModel<T>>(
    m: &M,
    spot: T,
    strike: T,
    rate: T,
    div: T,
    expiry: T,
    kind: OptionType,
    n_terms: usize,
) -> T {
    if expiry.val() <= 0.0 {
        return if kind.is_call() {
            (spot - strike).max_val(T::zero())
        } else {
            (strike - spot).max_val(T::zero())
        };
    }
    let (c1, c2, c4) = m.cumulants(rate, div, expiry);
    let width = T::from_f64(12.0) * (c2 + c4.max_val(T::zero()).sqrt()).sqrt();
    if width.val() < 1e-12 {
        let forward = spot * ((rate - div) * expiry).exp();
        let intrinsic = if kind.is_call() {
            (forward - strike).max_val(T::zero())
        } else {
            (strike - forward).max_val(T::zero())
        };
        return (-rate * expiry).exp() * intrinsic;
    }
    let a = c1 - width;
    let b = c1 + width;
    let w = b - a;
    let mut sum = T::zero();
    for n in 0..n_terms {
        let u = T::from_f64(n as f64 * std::f64::consts::PI) / w;
        let (re, im) = m.char_exponent(u);
        let drift = rate - div + m.martingale_compensator();
        let phire = (re * expiry).exp() * (im * expiry + u * drift * expiry).cos();
        let phiim = (re * expiry).exp() * (im * expiry + u * drift * expiry).sin();
        let x0 = (spot / strike).ln();
        let phase = u * (x0 - a);
        let real = phire * phase.cos() - phiim * phase.sin();
        let term = real * coeff_put(a, w, n);
        sum += if n == 0 { term * 0.5 } else { term };
    }
    let put = strike * (-rate * expiry).exp() * sum;
    if kind.is_call() {
        put + spot * (-div * expiry).exp() - strike * (-rate * expiry).exp()
    } else {
        put
    }
}
#[allow(clippy::too_many_arguments)]
pub fn cos_levy_price<T: Real, M: LevyModel<T>>(
    model: &M,
    spot: T,
    strike: T,
    rate: T,
    div_yield: T,
    expiry: T,
    option_type: OptionType,
    n_terms: usize,
) -> T {
    price_one(
        model,
        spot,
        strike,
        rate,
        div_yield,
        expiry,
        option_type,
        n_terms,
    )
}
#[allow(clippy::too_many_arguments)]
pub fn cos_levy_price_strip<T: Real, M: LevyModel<T>>(
    model: &M,
    spot: T,
    strikes: &[T],
    rate: T,
    div_yield: T,
    expiry: T,
    option_type: OptionType,
    n_terms: usize,
) -> Vec<T> {
    if expiry.val() <= 0.0 {
        return strikes
            .iter()
            .map(|&k| {
                if option_type.is_call() {
                    (spot - k).max_val(T::zero())
                } else {
                    (k - spot).max_val(T::zero())
                }
            })
            .collect();
    }
    let (c1, c2, c4) = model.cumulants(rate, div_yield, expiry);
    let width = T::from_f64(12.0) * (c2 + c4.max_val(T::zero()).sqrt()).sqrt();
    if width.val() < 1e-12 {
        let forward = spot * ((rate - div_yield) * expiry).exp();
        let discount = (-rate * expiry).exp();
        return strikes
            .iter()
            .map(|&k| {
                let intrinsic = if option_type.is_call() {
                    (forward - k).max_val(T::zero())
                } else {
                    (k - forward).max_val(T::zero())
                };
                discount * intrinsic
            })
            .collect();
    }
    let a = c1 - width;
    let w = width * 2.0;
    let drift = rate - div_yield + model.martingale_compensator();
    let mut phi = Vec::with_capacity(n_terms);
    for n in 0..n_terms {
        let u = T::from_f64(n as f64 * std::f64::consts::PI) / w;
        let (re, im) = model.char_exponent(u);
        let amp = (re * expiry).exp();
        let phase = im * expiry + u * drift * expiry;
        phi.push((amp * phase.cos(), amp * phase.sin()));
    }
    let coeffs: Vec<T> = (0..n_terms).map(|n| coeff_put(a, w, n)).collect();
    strikes
        .iter()
        .map(|&k| {
            let x0 = (spot / k).ln();
            let mut sum = T::zero();
            for n in 0..n_terms {
                let u = T::from_f64(n as f64 * std::f64::consts::PI) / w;
                let phase = u * (x0 - a);
                let (pr, pi) = phi[n];
                let term = (pr * phase.cos() - pi * phase.sin()) * coeffs[n];
                sum += if n == 0 { term * 0.5 } else { term };
            }
            let put = k * (-rate * expiry).exp() * sum;
            if option_type.is_call() {
                put + spot * (-div_yield * expiry).exp() - k * (-rate * expiry).exp()
            } else {
                put
            }
        })
        .collect()
}
