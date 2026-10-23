use ql_ad::{Real, ift_lift_1d};
#[inline]
pub fn erf<T: Real>(x: T) -> T {
    x.erf()
}
#[inline]
pub fn erfc<T: Real>(x: T) -> T {
    x.erfc()
}
#[inline]
pub fn norm_cdf<T: Real>(x: T) -> T {
    x.norm_cdf()
}
#[inline]
pub fn norm_pdf<T: Real>(x: T) -> T {
    x.norm_pdf()
}

/// Standard bivariate normal density with correlation `rho`.
pub fn bivariate_normal_pdf<T: Real>(a: T, b: T, rho: T) -> T {
    let one_minus = (T::one() - rho * rho).max_val(T::from_f64(1e-15));
    let exponent =
        -(a * a - T::from_f64(2.0) * rho * a * b + b * b) / (T::from_f64(2.0) * one_minus);
    exponent.exp() / (T::from_f64(2.0 * std::f64::consts::PI) * one_minus.sqrt())
}

/// Bivariate standard-normal CDF using the 20-point Genz/Drezner rule.
/// The quadrature nodes and weights are constants; all arithmetic involving
/// inputs remains in `T`, preserving Dual/Tape derivatives.
pub fn bivariate_normal_cdf<T: Real>(a: T, b: T, rho: T) -> T {
    let pa = a.norm_cdf();
    let pb = b.norm_cdf();
    if rho.val() >= 1.0 - 1e-12 {
        return pa.min_val(pb);
    }
    if rho.val() <= -1.0 + 1e-12 {
        return (pa + pb - T::one()).max_val(T::zero());
    }
    if rho.val().abs() < 1e-15 {
        return pa * pb;
    }
    const X: [f64; 10] = [
        0.9931285991850949,
        0.9639719272779138,
        0.912_234_428_251_326,
        0.8391169718222188,
        0.7463319064601508,
        0.636_053_680_726_515,
        0.5108670019508271,
        0.37370608871541955,
        0.2277858511416451,
        0.07652652113349733,
    ];
    const W: [f64; 10] = [
        0.017614007139152118,
        0.04060142980038694,
        0.06267204833410906,
        0.08327674157670475,
        0.1019301198172404,
        0.1181945319615184,
        0.13168863844917663,
        0.14209610931838205,
        0.14917298647260375,
        0.15275338713072585,
    ];
    if rho.val().abs() >= 0.925 {
        let h = -a;
        let mut k = -b;
        let mut hk = h * k;
        let negative = rho.val() < 0.0;
        if negative {
            k = -k;
            hk = -hk;
        }
        let ass = (T::one() - rho) * (T::one() + rho);
        let aa = ass.sqrt();
        let (bsqrt, max_hk) = if h.val() >= k.val() {
            (h - k, h)
        } else {
            (k - h, k)
        };
        let bs = bsqrt * bsqrt;
        let c = (T::from_f64(4.0) - hk) / 8.0;
        let d = (T::from_f64(12.0) - hk) / 16.0;
        let asr = -(bs / ass + hk) * 0.5;
        let mut value = T::zero();
        if asr.val() > -100.0 {
            value += aa
                * asr.exp()
                * (T::one() - c * (bs - ass) * (T::one() - d * bs / 5.0) / 3.0
                    + c * d * ass * ass / 5.0);
        }
        if (-hk).val() < 100.0 {
            value -= (-hk * 0.5).exp()
                * T::from_f64((2.0 * std::f64::consts::PI).sqrt())
                * (-bsqrt / aa).norm_cdf()
                * bsqrt
                * (T::one() - c * bs * (T::one() - d * bs / 5.0) / 3.0);
        }
        let half_a = aa * 0.5;
        for i in 0..10 {
            for &sign in &[-1.0_f64, 1.0_f64] {
                let xs = half_a * (1.0 - sign * X[i]);
                let xs2 = xs * xs;
                let rs = (T::one() - xs2).sqrt();
                let ar = -(bs / xs2 + hk) * 0.5;
                if ar.val() > -100.0 {
                    value += T::from_f64(W[i])
                        * half_a
                        * ar.exp()
                        * (((-hk * (T::one() - rs) / (T::from_f64(2.0) * (T::one() + rs))).exp()
                            / rs)
                            - (T::one() + c * xs2 * (T::one() + d * xs2)));
                }
            }
        }
        value /= -T::from_f64(2.0 * std::f64::consts::PI);
        if !negative {
            value += (-max_hk).norm_cdf();
        } else {
            value = -value;
            if k.val() > h.val() {
                value += if h.val() >= 0.0 {
                    (-h).norm_cdf() - (-k).norm_cdf()
                } else {
                    k.norm_cdf() - h.norm_cdf()
                };
            }
        }
        return value;
    }
    // asin(rho) expressed through Real's atan/sqrt primitives.
    let angle = (rho / (T::one() - rho * rho).sqrt()).atan();
    let mut sum = T::zero();
    for i in 0..10 {
        for &sign in &[-1.0_f64, 1.0_f64] {
            let x = sign * X[i];
            let theta = angle * (T::one() - x) * 0.5;
            let s = theta.sin();
            let denom = T::one() - s * s;
            let exponent = (s * a * b - T::from_f64(0.5) * (a * a + b * b)) / denom;
            sum += exponent.exp() * W[i];
        }
    }
    pa * pb + angle * sum / T::from_f64(4.0 * std::f64::consts::PI)
}
pub fn inv_norm_cdf_f64(p: f64) -> f64 {
    assert!(p > 0.0 && p < 1.0);
    let a = [
        -3.969683028665376e1,
        2.209460984245205e2,
        -2.759285104469687e2,
        1.38357751867269e2,
        -3.066479806614716e1,
        2.506628277459239,
    ];
    let b = [
        -5.447609879822406e1,
        1.615858368580409e2,
        -1.556989798598866e2,
        6.680131188771972e1,
        -1.328068155288572e1,
    ];
    let c = [
        -7.784894002430293e-3,
        -3.223964580411365e-1,
        -2.400758277161838,
        -2.549732539343734,
        4.374664141464968,
        2.938163982698783,
    ];
    let d = [
        7.784695709041462e-3,
        3.224671290700398e-1,
        2.445134137142996,
        3.754408661907416,
    ];
    let pl = 0.02425;
    let ph = 1.0 - pl;
    let mut x;
    if p < pl {
        let q = (-2.0 * p.ln()).sqrt();
        let num = ((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5];
        let den = (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0;
        x = num / den;
    } else if p > ph {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        let num = ((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5];
        let den = (((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0;
        x = -num / den;
    } else {
        let q = p - 0.5;
        let r = q * q;
        let num = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q;
        let den = ((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0;
        x = num / den;
    }
    for _ in 0..2 {
        let e = norm_cdf_f64(x) - p;
        let pdf = (-0.5 * x * x).exp() / (2.0 * std::f64::consts::PI).sqrt();
        if pdf > 0.0 {
            x -= e / pdf;
        }
    }
    x
}

pub fn inv_norm_cdf<T: Real>(p: T) -> T {
    let root = inv_norm_cdf_f64(p.val());
    let pdf = (-0.5 * root * root).exp() / (2.0 * std::f64::consts::PI).sqrt();
    let fxx = -root * pdf;
    ift_lift_1d(root, |x: T| x.norm_cdf() - p, T::from_f64(pdf), fxx)
}

fn norm_cdf_f64(x: f64) -> f64 {
    0.5 * libm::erfc(-x * std::f64::consts::FRAC_1_SQRT_2)
}
