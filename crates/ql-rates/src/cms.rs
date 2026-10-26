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

#[derive(Clone, Debug)]
pub struct MarkovFunctionalSwaptionSpec<T: Real = f64> {
    pub expiry: f64,
    pub payment_times: Vec<f64>,
    pub mean_reversion: T,
    pub atm_normal_vol: T,
    pub sabr_rho: T,
    pub sabr_nu: T,
    pub strike: T,
    pub notional: T,
    pub is_payer: bool,
}

#[derive(Copy, Clone, Debug)]
pub struct MarkovFunctionalValuation<T: Real> {
    pub swaption_npv: T,
    pub digital_call_pv: T,
    pub cms_swaplet_rate: T,
    pub forward_swap_rate: T,
    pub annuity: T,
}

pub fn price_markov_functional_1f<T: Real>(
    spec: &MarkovFunctionalSwaptionSpec<T>,
    curve: &impl YieldTermStructure<T>,
    n_gauss_hermite: usize,
) -> MarkovFunctionalValuation<T> {
    let te = spec.expiry.max(1e-6);
    let mut ann0 = T::zero();
    let mut prev_t = te;
    let mut last_t = te;
    for &t in &spec.payment_times {
        if t <= te {
            continue;
        }
        let tau = (t - prev_t).max(0.0);
        ann0 += curve.discount_factor(t) * tau;
        prev_t = t;
        last_t = t;
    }
    if ann0.val() <= 1e-14 || last_t <= te {
        return MarkovFunctionalValuation {
            swaption_npv: T::zero(),
            digital_call_pv: T::zero(),
            cms_swaplet_rate: T::zero(),
            forward_swap_rate: T::zero(),
            annuity: T::zero(),
        };
    }
    let p_te = curve.discount_factor(te);
    let p_tn = curve.discount_factor(last_t);
    let s0 = (p_te - p_tn) / ann0;

    let n_q = n_gauss_hermite.clamp(8, 32);
    let (gl_x, gl_w) = gauss_legendre(n_q);
    let mut z_nodes = Vec::with_capacity(n_q);
    let mut n_weights = Vec::with_capacity(n_q);
    let mut w_norm = 0.0;
    for (&x, &w) in gl_x.iter().zip(gl_w.iter()) {
        let z = 4.5 * x;
        let pdf = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
        let nw = 4.5 * w * pdf;
        z_nodes.push(z);
        n_weights.push(nw);
        w_norm += nw;
    }
    for nw in &mut n_weights {
        *nw /= w_norm.max(1e-14);
    }

    let g_hw = |t1: f64, t2: f64| -> T {
        let tau = (t2 - t1).max(0.0);
        let kt = spec.mean_reversion * tau;
        if kt.val().abs() < 1e-5 {
            T::from_f64(tau) - spec.mean_reversion * (0.5 * tau * tau)
        } else {
            (T::one() - (-kt).exp()) / spec.mean_reversion
        }
    };
    let g_tn = g_hw(te, last_t);
    let sig_x = spec.atm_normal_vol * te.sqrt();

    let mut defl_ann = vec![T::zero(); n_q];
    let mut prev = te;
    for &t in &spec.payment_times {
        if t <= te {
            continue;
        }
        let tau = (t - prev).max(0.0);
        let p0_ratio = curve.discount_factor(t) / p_tn.max_val(T::from_f64(1e-12));
        let dg = g_hw(te, t) - g_tn;
        let mut raw_bond = vec![T::zero(); n_q];
        let mut mean_bond = T::zero();
        for q in 0..n_q {
            let b_q = (-dg * sig_x * z_nodes[q]).exp();
            raw_bond[q] = b_q;
            mean_bond += b_q * n_weights[q];
        }
        let scale = p0_ratio / mean_bond.max_val(T::from_f64(1e-12));
        for q in 0..n_q {
            defl_ann[q] += raw_bond[q] * scale * tau;
        }
        prev = t;
    }

    let mut sum_ann_w = T::zero();
    for q in 0..n_q {
        sum_ann_w += defl_ann[q] * n_weights[q];
    }
    let mut w_ann = vec![T::zero(); n_q];
    for q in 0..n_q {
        w_ann[q] = defl_ann[q] * n_weights[q] / sum_ann_w.max_val(T::from_f64(1e-12));
    }

    let mut s_raw = vec![T::zero(); n_q];
    let mut mean_s_ann = T::zero();
    for q in 0..n_q {
        let z = z_nodes[q];
        let sr = s0
            + sig_x * z
            + spec.sabr_rho * spec.sabr_nu * sig_x * te.sqrt() * (0.5 * (z * z - 1.0))
            + spec.sabr_nu * spec.sabr_nu * sig_x * te * ((z * z * z - 3.0 * z) / 6.0);
        s_raw[q] = sr;
        mean_s_ann += w_ann[q] * sr;
    }
    let shift = s0 - mean_s_ann;

    let mut raw_cms_w = vec![T::zero(); n_q];
    let mut sum_cms_w = T::zero();
    for q in 0..n_q {
        let cw = (g_tn * sig_x * z_nodes[q]).exp() * n_weights[q];
        raw_cms_w[q] = cw;
        sum_cms_w += cw;
    }

    let eps = T::from_f64(1e-6);
    let mut unit_swaption = T::zero();
    let mut unit_digital = T::zero();
    let mut cms_rate = T::zero();
    for q in 0..n_q {
        let s_q = s_raw[q] + shift;
        let diff = if spec.is_payer {
            s_q - spec.strike
        } else {
            spec.strike - s_q
        };
        let payoff = (diff + (diff * diff + eps * eps).sqrt()) * 0.5;
        unit_swaption += w_ann[q] * payoff;
        let dig_diff = s_q - spec.strike;
        let dig =
            T::from_f64(0.5) * (T::one() + dig_diff / (dig_diff * dig_diff + eps * eps).sqrt());
        unit_digital += w_ann[q] * dig;
        let w_cms_q = raw_cms_w[q] / sum_cms_w.max_val(T::from_f64(1e-12));
        cms_rate += w_cms_q * s_q;
    }

    MarkovFunctionalValuation {
        swaption_npv: spec.notional * ann0 * unit_swaption,
        digital_call_pv: spec.notional * ann0 * unit_digital,
        cms_swaplet_rate: cms_rate,
        forward_swap_rate: s0,
        annuity: ann0,
    }
}

#[cfg(test)]
mod domain_suite_29_tests {
    use super::*;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;

    #[test]
    fn markov_functional_1f_parity_and_ad() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let c = FlatForwardCurve::new(d, 0.03);
        let mut s = MarkovFunctionalSwaptionSpec {
            expiry: 1.0,
            payment_times: vec![2.0, 3.0, 4.0, 5.0],
            mean_reversion: 0.05,
            atm_normal_vol: 0.008,
            sabr_rho: -0.25,
            sabr_nu: 0.35,
            strike: 0.031,
            notional: 100.0,
            is_payer: true,
        };
        let pay = price_markov_functional_1f(&s, &c, 16);
        s.is_payer = false;
        let rec = price_markov_functional_1f(&s, &c, 16);
        let intrinsic = s.notional * pay.annuity * (pay.forward_swap_rate - s.strike);
        assert!((pay.swaption_npv - rec.swaption_npv - intrinsic).abs() < 1e-11);
        assert!(pay.cms_swaplet_rate > pay.forward_swap_rate);

        let h = 1e-5;
        let eval_v = |v: f64| {
            let mut sp = s.clone();
            sp.is_payer = true;
            sp.atm_normal_vol = v;
            price_markov_functional_1f(&sp, &c, 16).swaption_npv
        };
        let fd = (eval_v(0.008 + h) - eval_v(0.008 - h)) / (2.0 * h);
        type D = Dual<1>;
        let cd = FlatForwardCurve::new(d, D::from_value(0.03));
        let sd = MarkovFunctionalSwaptionSpec {
            expiry: 1.0,
            payment_times: vec![2.0, 3.0, 4.0, 5.0],
            mean_reversion: D::from_value(0.05),
            atm_normal_vol: D::variable(0.008, 0),
            sabr_rho: D::from_value(-0.25),
            sabr_nu: D::from_value(0.35),
            strike: D::from_value(0.031),
            notional: D::from_value(100.0),
            is_payer: true,
        };
        let ad = price_markov_functional_1f(&sd, &cd, 16).swaption_npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
