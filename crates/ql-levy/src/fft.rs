use crate::LevyModel;
use ql_ad::Real;
use ql_instruments::OptionType;

#[allow(clippy::too_many_arguments)]
pub fn lewis_levy_price_strip<T: Real, M: LevyModel<T>>(
    model: &M,
    spot: T,
    strikes: &[T],
    rate: T,
    div_yield: T,
    expiry: T,
    option_type: OptionType,
    n_quad: usize,
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
    let n = n_quad.max(128);
    let u_max = 120.0_f64;
    let du = u_max / n as f64;
    let drift = rate - div_yield + model.martingale_compensator();
    let beta = T::from_f64(0.5);
    let mut phi = Vec::with_capacity(n);
    for j in 0..n {
        let u_f = (j as f64 + 0.5) * du;
        let u = T::from_f64(u_f);
        let (re_psi, im_psi) = model.char_exponent_shifted(u, beta);
        let amp = ((drift * 0.5 + re_psi) * expiry).exp() * T::from_f64(du / (u_f * u_f + 0.25));
        let ph = (u * drift + im_psi) * expiry;
        phi.push((u, amp * ph.cos(), amp * ph.sin()));
    }
    let df_r = (-rate * expiry).exp();
    let df_q = (-div_yield * expiry).exp();
    let inv_pi = T::from_f64(1.0 / std::f64::consts::PI);
    strikes
        .iter()
        .map(|&k| {
            let x = (spot / k).ln();
            let mut integral = T::zero();
            for &(u, pr, pi) in &phi {
                let ux = u * x;
                integral += pr * ux.cos() - pi * ux.sin();
            }
            let call =
                (spot * df_q - (spot * k).sqrt() * df_r * inv_pi * integral).max_val(T::zero());
            if option_type.is_call() {
                call
            } else {
                (call - spot * df_q + k * df_r).max_val(T::zero())
            }
        })
        .collect()
}

fn radix2_fft_in_place<T: Real>(buf: &mut [(T, T)]) {
    let n = buf.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            buf.swap(i, j);
        }
    }
    let mut len = 2usize;
    while len <= n {
        let half = len / 2;
        let ang = -std::f64::consts::TAU / len as f64;
        for start in (0..n).step_by(len) {
            for k in 0..half {
                let theta = ang * k as f64;
                let wr = T::from_f64(theta.cos());
                let wi = T::from_f64(theta.sin());
                let (ur, ui) = buf[start + k];
                let (vr, vi) = buf[start + k + half];
                let tr = vr * wr - vi * wi;
                let ti = vr * wi + vi * wr;
                buf[start + k] = (ur + tr, ui + ti);
                buf[start + k + half] = (ur - tr, ui - ti);
            }
        }
        len <<= 1;
    }
}

#[allow(clippy::too_many_arguments)]
pub fn carr_madan_fft_levy_price_strip<T: Real, M: LevyModel<T>>(
    model: &M,
    spot: T,
    strikes: &[T],
    rate: T,
    div_yield: T,
    expiry: T,
    option_type: OptionType,
    fft_power_of_two: usize,
    damping_alpha: f64,
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
    let p = fft_power_of_two.clamp(8, 13);
    let n = 1usize << p;
    let alpha = damping_alpha.clamp(0.25, 4.0);
    let eta = 0.15_f64;
    let lambda = std::f64::consts::TAU / (n as f64 * eta);
    let b = 0.5 * n as f64 * lambda;
    let drift = rate - div_yield + model.martingale_compensator();
    let beta = T::from_f64(alpha + 1.0);
    let df_r = (-rate * expiry).exp();
    let df_q = (-div_yield * expiry).exp();
    let mut buf = vec![(T::zero(), T::zero()); n];
    for (j, slot) in buf.iter_mut().enumerate() {
        let v_f = j as f64 * eta;
        let v = T::from_f64(v_f);
        let (re_psi, im_psi) = model.char_exponent_shifted(v, beta);
        let amp = df_r * ((drift * beta + re_psi) * expiry).exp();
        let ph = (v * drift + im_psi) * expiry + T::from_f64(b * v_f);
        let phi_r = amp * ph.cos();
        let phi_i = amp * ph.sin();
        let den_r = alpha * alpha + alpha - v_f * v_f;
        let den_i = (2.0 * alpha + 1.0) * v_f;
        let den_norm = den_r * den_r + den_i * den_i;
        let simpson = eta / 3.0
            * if j == 0 {
                1.0
            } else if j % 2 == 1 {
                4.0
            } else {
                2.0
            };
        let scale = T::from_f64(simpson / den_norm);
        let qr = (phi_r * T::from_f64(den_r) + phi_i * T::from_f64(den_i)) * scale;
        let qi = (phi_i * T::from_f64(den_r) - phi_r * T::from_f64(den_i)) * scale;
        *slot = (qr, qi);
    }
    radix2_fft_in_place(&mut buf);
    let inv_pi = 1.0 / std::f64::consts::PI;
    let grid_calls: Vec<T> = (0..n)
        .map(|m| {
            let km = -b + m as f64 * lambda;
            buf[m].0 * T::from_f64((-alpha * km).exp() * inv_pi)
        })
        .collect();
    strikes
        .iter()
        .map(|&k| {
            let log_k = (k / spot).ln();
            let u = (log_k + T::from_f64(b)) / T::from_f64(lambda);
            let idx = (u.val().floor() as isize).clamp(1, n as isize - 3) as usize;
            let t = u - T::from_f64(idx as f64);
            let y0 = grid_calls[idx - 1];
            let y1 = grid_calls[idx];
            let y2 = grid_calls[idx + 1];
            let y3 = grid_calls[idx + 2];
            let c0 = y1;
            let c1 = (y2 - y0) * 0.5;
            let c2 = y0 - y1 * 2.5 + y2 * 2.0 - y3 * 0.5;
            let c3 = (y3 - y0) * 0.5 + (y1 - y2) * 1.5;
            let call = (spot * (c0 + t * (c1 + t * (c2 + t * c3)))).max_val(T::zero());
            if option_type.is_call() {
                call
            } else {
                (call - spot * df_q + k * df_r).max_val(T::zero())
            }
        })
        .collect()
}

#[cfg(test)]
mod domain_suite_26_tests {
    use super::*;
    use crate::{VarianceGamma, cos_levy_price_strip};
    use ql_ad::Dual;

    #[test]
    fn lewis_and_carr_madan_fft_match_cos_strip_and_ad() {
        let vg = VarianceGamma {
            sigma: 0.20,
            nu: 0.25,
            theta: -0.10,
        };
        let strikes = vec![90.0, 95.0, 100.0, 105.0, 110.0];
        let cos =
            cos_levy_price_strip(&vg, 100.0, &strikes, 0.03, 0.01, 1.0, OptionType::Call, 256);
        let lewis =
            lewis_levy_price_strip(&vg, 100.0, &strikes, 0.03, 0.01, 1.0, OptionType::Call, 512);
        let fft = carr_madan_fft_levy_price_strip(
            &vg,
            100.0,
            &strikes,
            0.03,
            0.01,
            1.0,
            OptionType::Call,
            11,
            1.5,
        );
        for i in 0..strikes.len() {
            assert!(
                (lewis[i] - cos[i]).abs() < 1e-3,
                "lewis vs cos mismatch at {}: {} vs {}",
                strikes[i],
                lewis[i],
                cos[i]
            );
            assert!(
                (fft[i] - cos[i]).abs() < 2e-3,
                "fft vs cos mismatch at {}: {} vs {}",
                strikes[i],
                fft[i],
                cos[i]
            );
        }
        type D = Dual<1>;
        let vg_d = VarianceGamma {
            sigma: D::variable(0.20, 0),
            nu: D::from_value(0.25),
            theta: D::from_value(-0.10),
        };
        let ks_d = vec![D::from_value(100.0)];
        let d_fft = carr_madan_fft_levy_price_strip(
            &vg_d,
            D::from_value(100.0),
            &ks_d,
            D::from_value(0.03),
            D::from_value(0.01),
            D::from_value(1.0),
            OptionType::Call,
            10,
            1.5,
        )[0];
        assert!(d_fft.d[0] > 10.0);
    }
}
