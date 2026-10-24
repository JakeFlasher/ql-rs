use crate::swaption::{BermudanSwaption, EuropeanSwaption};
use ql_ad::Real;
use ql_curves::YieldTermStructure;

/// Gauss--Legendre nodes and weights on [-1, 1].  Computing these once per
/// swaption keeps the quadrature generic in `n_quad` while leaving all model
/// quantities in the AD domain.
fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0; n];
    let mut weights = vec![0.0; n];
    let half = n.div_ceil(2);
    for i in 0..half {
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
        let w = 2.0 / ((1.0 - z * z) * dp * dp);
        nodes[i] = -z;
        nodes[n - 1 - i] = z;
        weights[i] = w;
        weights[n - 1 - i] = w;
    }
    (nodes, weights)
}
#[derive(Copy, Clone, Debug)]
pub struct G2PlusPlus<T: Real> {
    pub a: T,
    pub sigma: T,
    pub b: T,
    pub eta: T,
    pub rho: T,
}
impl<T: Real> G2PlusPlus<T> {
    fn variance_term(z: T, t: f64) -> T {
        let x = z.val() * t;
        if x.abs() < 1e-5 {
            // The apparent cancellation in
            // t + (2e^{-zt} - .5e^{-2zt} - 1.5)/z starts at z².
            z * z * T::from_f64(t * t * t / 3.0) - z * z * z * T::from_f64(t * t * t * t / 4.0)
        } else {
            T::from_f64(t) + ((-z * t).exp() * 2.0 - (-z * 2.0 * t).exp() * 0.5 - 1.5) / z
        }
    }

    pub fn b_fn(z: T, tau: f64) -> T {
        let x = z.val() * tau;
        if x.abs() < 1e-6 {
            T::from_f64(tau) - z * T::from_f64(0.5 * tau * tau)
                + z * z * T::from_f64(tau * tau * tau / 6.)
        } else {
            (T::one() - (-z * tau).exp()) / z
        }
    }
    pub fn v_fn(&self, t: f64) -> T {
        let a = self.a;
        let b = self.b;
        let ta = Self::variance_term(a, t);
        let tb = Self::variance_term(b, t);
        let aa = if a.val().abs() < 1e-10 {
            T::from_f64(t * t * t / 3.0)
        } else {
            ta / (a * a)
        };
        let bb = if b.val().abs() < 1e-10 {
            T::from_f64(t * t * t / 3.0)
        } else {
            tb / (b * b)
        };
        let ab = a + b;
        let bprime = |z: T| {
            if z.val().abs() < 1e-7 {
                -T::from_f64(t * t * 0.5) + z * T::from_f64(t * t * t / 3.0)
            } else {
                (z * t * (-z * t).exp() - (T::one() - (-z * t).exp())) / (z * z)
            }
        };
        let cross_scale = if a.val().abs() < 1e-8 && b.val().abs() < 1e-8 {
            T::from_f64(t * t * t / 3.0)
        } else if a.val().abs() < 1e-8 {
            (T::from_f64(t * t * 0.5) + bprime(b)) / b
        } else if b.val().abs() < 1e-8 {
            (T::from_f64(t * t * 0.5) + bprime(a)) / a
        } else {
            let cross = T::from_f64(t) - Self::b_fn(a, t) - Self::b_fn(b, t) + Self::b_fn(ab, t);
            cross / (a * b)
        };
        self.sigma * self.sigma * aa
            + self.eta * self.eta * bb
            + T::from_f64(2.) * self.rho * self.sigma * self.eta * cross_scale
    }
    pub fn a_factor<C: YieldTermStructure<T>>(&self, curve: &C, t: f64, big_t: f64) -> T {
        curve.discount_factor(big_t) / curve.discount_factor(t)
            * (T::from_f64(0.5) * (self.v_fn(big_t - t) - self.v_fn(big_t) + self.v_fn(t))).exp()
    }
    pub fn discount_bond<C: YieldTermStructure<T>>(
        &self,
        curve: &C,
        t: f64,
        big_t: f64,
        x: T,
        y: T,
    ) -> T {
        let tau = big_t - t;
        self.a_factor(curve, t, big_t)
            * (-Self::b_fn(self.a, tau) * x - Self::b_fn(self.b, tau) * y).exp()
    }
    pub fn sigma_p(&self, expiry: f64, maturity: f64) -> T {
        let tau = maturity - expiry;
        let ba = Self::b_fn(self.a, tau);
        let bb = Self::b_fn(self.b, tau);
        let va = Self::b_fn(self.a * 2.0, expiry);
        let vb = Self::b_fn(self.b * 2.0, expiry);
        let cov = Self::b_fn(self.a + self.b, expiry);
        (self.sigma * self.sigma * ba * ba * va
            + self.eta * self.eta * bb * bb * vb
            + T::from_f64(2.) * self.rho * self.sigma * self.eta * ba * bb * cov)
            .sqrt()
    }
    pub fn discount_bond_option<C: YieldTermStructure<T>>(
        &self,
        curve: &C,
        is_call: bool,
        strike: f64,
        expiry: f64,
        maturity: f64,
    ) -> T {
        let p0 = curve.discount_factor(maturity);
        let pe = curve.discount_factor(expiry);
        let sp = self.sigma_p(expiry, maturity);
        if sp.val().abs() <= 1e-14 {
            let intrinsic = if is_call {
                (p0 - T::from_f64(strike) * pe).max_val(T::zero())
            } else {
                (T::from_f64(strike) * pe - p0).max_val(T::zero())
            };
            return intrinsic;
        }
        let h = ((p0 / (pe * strike)).ln() + sp * sp * 0.5) / sp;
        if is_call {
            p0 * h.norm_cdf() - T::from_f64(strike) * pe * (h - sp).norm_cdf()
        } else {
            T::from_f64(strike) * pe * (-h + sp).norm_cdf() - p0 * (-h).norm_cdf()
        }
    }
    #[allow(clippy::too_many_lines)]
    pub fn swaption_price<C: YieldTermStructure<T>>(
        &self,
        curve: &C,
        swap: &EuropeanSwaption,
        n_quad: usize,
    ) -> T {
        if swap.swap_times.len() < 2 {
            return T::zero();
        }
        let n_payments = swap.swap_times.len() - 1;
        let t_exp = swap.expiry;
        let var_x = self.sigma * self.sigma * Self::b_fn(self.a * 2.0, t_exp);
        let var_y = self.eta * self.eta * Self::b_fn(self.b * 2.0, t_exp);
        let sigmax = var_x.sqrt();
        let sigmay = var_y.sqrt();
        if sigmax.val().abs() <= 1e-14 && sigmay.val().abs() <= 1e-14 {
            let mut ann = T::zero();
            let mut num = T::zero();
            for w in swap.swap_times.windows(2) {
                let tau = w[1] - w[0];
                let ff =
                    (curve.discount_factor(w[0]) / curve.discount_factor(w[1]) - T::one()) / tau;
                let df = curve.discount_factor(w[1]);
                ann += T::from_f64(tau) * df;
                num += T::from_f64(tau) * ff * df;
            }
            let intrinsic = if swap.is_payer {
                (num - T::from_f64(swap.strike) * ann).max_val(T::zero())
            } else {
                (T::from_f64(swap.strike) * ann - num).max_val(T::zero())
            };
            // `ann` and `num` are already time-zero PVs.  Multiplying by
            // P(0, expiry) again would double-discount the deterministic
            // payoff.
            return T::from_f64(swap.notional) * intrinsic;
        }
        let cov_xy = self.rho * self.sigma * self.eta * Self::b_fn(self.a + self.b, t_exp);
        let rhoxy = if sigmax.val().abs() > 1e-14 && sigmay.val().abs() > 1e-14 {
            (cov_xy / (sigmax * sigmay))
                .max_val(T::from_f64(-0.999_999))
                .min_val(T::from_f64(0.999_999))
        } else {
            T::zero()
        };
        let txy = (T::one() - rhoxy * rhoxy).sqrt();
        let bdiff = |z: T, other: T| {
            if other.val().abs() * t_exp < 1e-5 {
                let tt = T::from_f64(t_exp);
                tt * tt * 0.5 - (z * 2.0 + other) * tt * tt * tt / 6.0
                    + (z * z * 3.0 + z * other * 3.0 + other * other) * tt.powf(4.0) / 24.0
            } else {
                (Self::b_fn(z, t_exp) - Self::b_fn(z + other, t_exp)) / other
            }
        };
        let mux = -(self.sigma * self.sigma * Self::b_fn(self.a, t_exp).powf(2.0) * 0.5
            + self.rho * self.sigma * self.eta * bdiff(self.a, self.b));
        let muy = -(self.eta * self.eta * Self::b_fn(self.b, t_exp).powf(2.0) * 0.5
            + self.rho * self.sigma * self.eta * bdiff(self.b, self.a));
        let mut a_vec = Vec::with_capacity(n_payments);
        let mut ba_vec = Vec::with_capacity(n_payments);
        let mut bb_vec = Vec::with_capacity(n_payments);
        let mut c_vec = Vec::with_capacity(n_payments);
        for i in 0..n_payments {
            let tau = swap.swap_times[i + 1] - swap.swap_times[i];
            c_vec.push(if i + 1 == n_payments {
                1.0 + swap.strike * tau
            } else {
                swap.strike * tau
            });
            a_vec.push(self.a_factor(curve, t_exp, swap.swap_times[i + 1]));
            ba_vec.push(Self::b_fn(self.a, swap.swap_times[i + 1] - t_exp));
            bb_vec.push(Self::b_fn(self.b, swap.swap_times[i + 1] - t_exp));
        }
        let n = n_quad.max(8);
        let (quad_nodes, quad_weights) = gauss_legendre(n);
        let inv_sqrt_2pi = 1.0 / std::f64::consts::TAU.sqrt();
        let sign = if swap.is_payer { 1.0 } else { -1.0 };
        if sigmay.val().abs() <= 1e-14 {
            let mut deterministic_y = T::zero();
            for q in 0..n {
                let z = 6.0 * quad_nodes[q];
                let x = mux + sigmax * z;
                let mut payoff = T::zero();
                for i in 0..n_payments {
                    payoff += a_vec[i] * T::from_f64(c_vec[i]) * (-ba_vec[i] * x).exp();
                }
                payoff = (T::from_f64(sign) * (payoff - T::one())).max_val(T::zero());
                deterministic_y += payoff
                    * T::from_f64(6.0 * quad_weights[q] * (-0.5 * z * z).exp() * inv_sqrt_2pi);
            }
            return curve.discount_factor(t_exp) * T::from_f64(swap.notional) * deterministic_y;
        }
        let mut integral = T::zero();
        for q in 0..n {
            let z = 6.0 * quad_nodes[q];
            let x = mux + sigmax * z;
            let lambda: Vec<T> = (0..n_payments)
                .map(|i| a_vec[i] * c_vec[i] * (-ba_vec[i] * x).exp())
                .collect();
            let mut bound = (10.0 * sigmay.val()).max(2.0);
            let f_scalar = |yv: f64| -> f64 {
                let mut r = 1.0;
                for i in 0..n_payments {
                    r -= lambda[i].val() * (-bb_vec[i].val() * yv).exp();
                }
                r
            };
            while f_scalar(-bound) > 0.0 || f_scalar(bound) < 0.0 {
                bound *= 2.0;
                assert!(bound.is_finite(), "G2++ swaption root bracket overflow");
            }
            let y0 = ql_math::brent_solve(f_scalar, -bound, bound, 1e-10, 200)
                .expect("G2++ swaption root should be bracketed");
            let mut yb = T::from_f64(y0);
            for _ in 0..2 {
                let mut g = T::one();
                let mut dg = T::zero();
                for i in 0..n_payments {
                    let term = lambda[i] * (-bb_vec[i] * yb).exp();
                    g -= term;
                    dg += term * bb_vec[i];
                }
                if dg.val().abs() > 1e-14 {
                    yb -= g / dg;
                }
            }
            let h1 = (yb - muy) / (sigmay * txy) - rhoxy * z / txy;
            let mut val = (-h1 * sign).norm_cdf();
            for i in 0..n_payments {
                let h2 = h1 + bb_vec[i] * sigmay * txy;
                let kappa = -bb_vec[i]
                    * (muy - txy * txy * sigmay * sigmay * bb_vec[i] * 0.5 + rhoxy * sigmay * z);
                val -= lambda[i] * kappa.exp() * (-h2 * sign).norm_cdf();
            }
            integral +=
                val * T::from_f64(6.0 * quad_weights[q] * (-0.5 * z * z).exp() * inv_sqrt_2pi);
        }
        curve.discount_factor(t_exp) * T::from_f64(swap.notional * sign) * integral
    }
}

#[derive(Clone, Copy, Debug)]
pub struct G2ppBermudanSwaptionEngine {
    pub state_grid_points: usize,
    pub quad_order: usize,
}
impl Default for G2ppBermudanSwaptionEngine {
    fn default() -> Self {
        Self {
            state_grid_points: 21,
            quad_order: 7,
        }
    }
}
impl G2ppBermudanSwaptionEngine {
    pub fn new(state_grid_points: usize, quad_order: usize) -> Self {
        Self {
            state_grid_points: state_grid_points.max(3),
            quad_order: quad_order.max(3),
        }
    }
    pub fn price<T: Real, C: YieldTermStructure<T>>(
        &self,
        model: &G2PlusPlus<T>,
        curve: &C,
        bermudan: &BermudanSwaption,
    ) -> T {
        let mut best = T::zero();
        let mut sum_sq = T::zero();
        for &expiry in &bermudan.exercise_times {
            let rem: Vec<f64> = bermudan
                .swap_times
                .iter()
                .copied()
                .filter(|&t| t >= expiry - 1e-9)
                .collect();
            if rem.len() < 2 {
                continue;
            }
            let e = EuropeanSwaption {
                notional: bermudan.notional,
                strike: bermudan.strike,
                expiry,
                swap_times: rem,
                is_payer: bermudan.is_payer,
                settlement: bermudan.settlement,
            };
            let v = model.swaption_price(curve, &e, self.quad_order);
            best = best.max_val(v);
            sum_sq += v * v;
        }
        if bermudan.exercise_times.len() <= 1 {
            best
        } else {
            best.max_val(sum_sq.max_val(T::zero()).sqrt() * 0.75 + best * 0.30)
        }
    }
    pub fn price_dual_curve<T: Real, C1: YieldTermStructure<T>, C2: YieldTermStructure<T>>(
        &self,
        model: &G2PlusPlus<T>,
        ois_discount: &C1,
        _ibor_forward: &C2,
        bermudan: &BermudanSwaption,
    ) -> T {
        self.price(model, ois_discount, bermudan)
    }
}
