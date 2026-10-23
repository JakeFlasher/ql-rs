use ql_ad::Real;
use ql_curves::YieldTermStructure;

fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0; n];
    let mut w = vec![0.0; n];
    for i in 0..n.div_ceil(2) {
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..32 {
            let mut p0 = 1.0;
            let mut p1 = z;
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            let dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
            let dz = p1 / dp;
            z -= dz;
            if dz.abs() < 1e-15 {
                break;
            }
        }
        let mut p0 = 1.0;
        let mut p1 = z;
        for k in 2..=n {
            let p2 = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
            p0 = p1;
            p1 = p2;
        }
        let dp = n as f64 * (z * p1 - p0) / (z * z - 1.0);
        let wi = 2.0 / ((1.0 - z * z) * dp * dp);
        x[i] = -z;
        x[n - 1 - i] = z;
        w[i] = wi;
        w[n - 1 - i] = wi;
    }
    (x, w)
}
#[derive(Copy, Clone, Debug)]
pub struct LinearTsrCmsPricer<T: Real> {
    pub mean_reversion: T,
}
impl<T: Real> LinearTsrCmsPricer<T> {
    fn g(&self, t: f64, tt: f64) -> T {
        let tau = (tt - t).max(0.);
        if (self.mean_reversion.val() * tau).abs() < 1e-6 {
            T::from_f64(tau) - self.mean_reversion * (tau * tau * 0.5)
                + self.mean_reversion * self.mean_reversion * (tau * tau * tau / 6.0)
        } else {
            (T::one() - (-self.mean_reversion * tau).exp()) / self.mean_reversion
        }
    }
    fn params<C: YieldTermStructure<T>>(
        &self,
        curve: &C,
        tfix: f64,
        tpay: f64,
        times: &[f64],
    ) -> (T, T, T, T) {
        let mut a = T::zero();
        let mut gx = T::zero();
        for w in times.windows(2) {
            let tau = w[1] - w[0];
            let p = curve.discount_factor(w[1]);
            a += T::from_f64(tau) * p;
            gx += T::from_f64(tau) * p * self.g(tfix, w[1]);
        }
        let s =
            (curve.discount_factor(times[0]) - curve.discount_factor(*times.last().unwrap())) / a;
        let gamma = gx / a;
        let aa = curve.discount_factor(tpay) * (gamma - self.g(tfix, tpay))
            / (curve.discount_factor(*times.last().unwrap())
                * self.g(tfix, *times.last().unwrap())
                + s * a * gamma);
        let bb = curve.discount_factor(tpay) / a - aa * s;
        (s, a, aa, bb)
    }
    fn psi(&self, s: T, k: T, tf: f64, vol: T, call: bool) -> T {
        let st = vol * tf.sqrt();
        let om = if call { 1. } else { -1. };
        let d = if st.val() > 1e-12 {
            (s - k) / st
        } else {
            T::zero()
        };
        T::from_f64(om) * (s - k) * (T::from_f64(om) * d).norm_cdf() + st * d.norm_pdf()
    }
    fn otm_psi(&self, s: T, k: T, tf: f64, vol: T) -> T {
        self.psi(s, k, tf, vol, k.val() >= s.val())
    }
    pub fn cms_caplet_price<C: YieldTermStructure<T>>(
        &self,
        curve: &C,
        tfix: f64,
        tpay: f64,
        times: &[f64],
        strike: T,
        smile: impl Fn(T) -> T,
    ) -> T {
        let (s, a, aa, bb) = self.params(curve, tfix, tpay, times);
        let vol = smile(s);
        let width = vol * tfix.sqrt() * 6.;
        let lo = strike;
        let hi = (s + width).max_val(strike);
        let mut integral = T::zero();
        let (nodes, weights) = gauss_legendre(32);
        for q in 0..32 {
            let z = 0.5 * (nodes[q] + 1.0);
            let k = lo + (hi - lo) * z;
            integral += self.otm_psi(s, k, tfix, smile(k)) * (hi - lo) * 0.5 * weights[q];
        }
        let singular = (s - strike).max_val(T::zero()) * (aa * s + bb)
            + (aa * strike + bb) * self.otm_psi(s, strike, tfix, vol);
        a * (singular + T::from_f64(2.) * aa * integral)
    }
    pub fn cms_floorlet_price<C: YieldTermStructure<T>>(
        &self,
        curve: &C,
        tfix: f64,
        tpay: f64,
        times: &[f64],
        strike: T,
        smile: impl Fn(T) -> T,
    ) -> T {
        let (s, a, aa, bb) = self.params(curve, tfix, tpay, times);
        let vol = smile(s);
        let width = vol * tfix.sqrt() * 6.;
        let lo = (s - width).min_val(strike);
        let hi = strike;
        let mut integral = T::zero();
        let (nodes, weights) = gauss_legendre(32);
        for q in 0..32 {
            let z = 0.5 * (nodes[q] + 1.0);
            let k = lo + (hi - lo) * z;
            integral += self.otm_psi(s, k, tfix, smile(k)) * (hi - lo) * 0.5 * weights[q];
        }
        let singular = (strike - s).max_val(T::zero()) * (aa * s + bb)
            + (aa * strike + bb) * self.otm_psi(s, strike, tfix, vol);
        a * (singular - T::from_f64(2.) * aa * integral)
    }
    pub fn cms_swaplet_rate<C: YieldTermStructure<T>>(
        &self,
        curve: &C,
        tfix: f64,
        tpay: f64,
        times: &[f64],
        smile: impl Fn(T) -> T,
    ) -> T {
        let (s, _, _, _) = self.params(curve, tfix, tpay, times);
        let cap = self.cms_caplet_price(curve, tfix, tpay, times, s, &smile);
        let floor = self.cms_floorlet_price(curve, tfix, tpay, times, s, smile);
        s + (cap - floor) / curve.discount_factor(tpay)
    }
}
