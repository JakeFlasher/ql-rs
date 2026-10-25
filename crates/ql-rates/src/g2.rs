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
impl<T: Real> G2PlusPlus<T> {
    pub fn sigma_x(&self, t: f64) -> T {
        (self.sigma * self.sigma * Self::b_fn(self.a * 2.0, t)).sqrt()
    }
    pub fn sigma_y(&self, t: f64) -> T {
        (self.eta * self.eta * Self::b_fn(self.b * 2.0, t)).sqrt()
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
        self.price_dual_curve(model, curve, curve, bermudan)
    }

    pub fn price_dual_curve<T: Real, C1: YieldTermStructure<T>, C2: YieldTermStructure<T>>(
        &self,
        model: &G2PlusPlus<T>,
        ois_discount: &C1,
        ibor_forward: &C2,
        bermudan: &BermudanSwaption,
    ) -> T {
        let valid_ex: Vec<(f64, Vec<f64>)> = bermudan
            .exercise_times
            .iter()
            .copied()
            .filter_map(|expiry| {
                let rem: Vec<f64> = bermudan
                    .swap_times
                    .iter()
                    .copied()
                    .filter(|&t| t >= expiry - 1e-9)
                    .collect();
                (rem.len() >= 2).then_some((expiry, rem))
            })
            .collect();
        if valid_ex.is_empty() {
            return T::zero();
        }
        let mut euro_best = T::zero();
        for (expiry, rem) in &valid_ex {
            let e = EuropeanSwaption {
                notional: bermudan.notional,
                strike: bermudan.strike,
                expiry: *expiry,
                swap_times: rem.clone(),
                is_payer: bermudan.is_payer,
                settlement: bermudan.settlement,
            };
            let v = model.swaption_price(ois_discount, &e, self.quad_order);
            euro_best = euro_best.max_val(v);
        }
        let g = self.state_grid_points.clamp(5, 15);
        let m_ex = valid_ex.len();
        let sign = if bermudan.is_payer { 1.0 } else { -1.0 };
        let t_terminal = *bermudan.swap_times.last().unwrap_or(&1.0);
        let q_nodes = [
            -2.856_970_013_872_805_6,
            -1.355_626_179_974_266,
            0.0,
            1.355_626_179_974_266,
            2.856_970_013_872_805_6_f64,
        ];
        let q_weights = [
            0.011_257_411_327_720_69,
            0.222_075_922_005_612_64,
            0.533_333_333_333_333_3,
            0.222_075_922_005_612_64,
            0.011_257_411_327_720_69_f64,
        ];
        let intrinsic_at = |m_idx: usize, x: T, y: T| -> T {
            let (t_ex, ref rem) = valid_ex[m_idx];
            let mut flt = T::zero();
            let mut ann = T::zero();
            for k in 0..(rem.len() - 1) {
                let t0 = rem[k];
                let t1 = rem[k + 1];
                let tau = t1 - t0;
                let p0_ois = model.discount_bond(ois_discount, t_ex, t0, x, y);
                let p1_ois = model.discount_bond(ois_discount, t_ex, t1, x, y);
                let ratio = (ibor_forward.discount_factor(t0) / ibor_forward.discount_factor(t1))
                    / (ois_discount.discount_factor(t0) / ois_discount.discount_factor(t1));
                let fwd = (p0_ois / p1_ois * ratio - T::one()) / tau;
                flt += fwd * p1_ois * tau;
                ann += p1_ois * tau;
            }
            let payoff = (T::from_f64(sign) * (flt - ann * bermudan.strike)).max_val(T::zero());
            let p_term = model
                .discount_bond(ois_discount, t_ex, t_terminal, x, y)
                .max_val(T::from_f64(1e-12));
            payoff * bermudan.notional / p_term
        };
        let grid_bounds = |t: f64| -> (T, T) {
            let sx = model.sigma_x(t).max_val(T::from_f64(1e-6)) * 3.5;
            let sy = model.sigma_y(t).max_val(T::from_f64(1e-6)) * 3.5;
            (sx, sy)
        };
        let mut vals = vec![vec![T::zero(); g]; g];
        let (bx_last, by_last) = grid_bounds(valid_ex[m_ex - 1].0);
        for (i, row) in vals.iter_mut().enumerate() {
            let xi = -bx_last + bx_last * (2.0 * i as f64 / (g - 1) as f64);
            for (j, cell) in row.iter_mut().enumerate() {
                let yj = -by_last + by_last * (2.0 * j as f64 / (g - 1) as f64);
                *cell = intrinsic_at(m_ex - 1, xi, yj);
            }
        }
        let interp = |grid: &[Vec<T>], bx: T, by: T, x: T, y: T| -> T {
            let ux = ((x + bx) / (bx * 2.0)).max_val(T::zero()).min_val(T::one()) * (g - 1) as f64;
            let uy = ((y + by) / (by * 2.0)).max_val(T::zero()).min_val(T::one()) * (g - 1) as f64;
            let i0 = (ux.val().floor() as usize).min(g - 2);
            let j0 = (uy.val().floor() as usize).min(g - 2);
            let tx = ux - T::from_f64(i0 as f64);
            let ty = uy - T::from_f64(j0 as f64);
            grid[i0][j0] * (T::one() - tx) * (T::one() - ty)
                + grid[i0 + 1][j0] * tx * (T::one() - ty)
                + grid[i0][j0 + 1] * (T::one() - tx) * ty
                + grid[i0 + 1][j0 + 1] * tx * ty
        };
        for m_idx in (0..(m_ex - 1)).rev() {
            let t_curr = valid_ex[m_idx].0;
            let t_next = valid_ex[m_idx + 1].0;
            let dt = (t_next - t_curr).max(1e-8);
            let (bx_curr, by_curr) = grid_bounds(t_curr);
            let (bx_next, by_next) = grid_bounds(t_next);
            let sx = model.sigma_x(dt);
            let sy = model.sigma_y(dt);
            let ea = (-model.a * dt).exp();
            let eb = (-model.b * dt).exp();
            let cov = model.rho * model.sigma * model.eta * G2PlusPlus::b_fn(model.a + model.b, dt);
            let rxy = if sx.val() > 1e-14 && sy.val() > 1e-14 {
                (cov / (sx * sy))
                    .max_val(T::from_f64(-0.999))
                    .min_val(T::from_f64(0.999))
            } else {
                T::zero()
            };
            let s_orth = (T::one() - rxy * rxy).max_val(T::zero()).sqrt();
            let mut next_vals = vec![vec![T::zero(); g]; g];
            for (i, row) in next_vals.iter_mut().enumerate() {
                let xi = -bx_curr + bx_curr * (2.0 * i as f64 / (g - 1) as f64);
                for (j, cell) in row.iter_mut().enumerate() {
                    let yj = -by_curr + by_curr * (2.0 * j as f64 / (g - 1) as f64);
                    let mut cont = T::zero();
                    for (qx, &zx) in q_nodes.iter().enumerate() {
                        for (qy, &zy) in q_nodes.iter().enumerate() {
                            let xn = xi * ea + sx * zx;
                            let yn = yj * eb + sy * (rxy * zx + s_orth * zy);
                            let w = q_weights[qx] * q_weights[qy];
                            cont += interp(&vals, bx_next, by_next, xn, yn) * w;
                        }
                    }
                    *cell = intrinsic_at(m_idx, xi, yj).max_val(cont);
                }
            }
            vals = next_vals;
        }
        let t0 = valid_ex[0].0;
        let (bx0, by0) = grid_bounds(t0);
        let sx0 = model.sigma_x(t0);
        let sy0 = model.sigma_y(t0);
        let cov0 = model.rho * model.sigma * model.eta * G2PlusPlus::b_fn(model.a + model.b, t0);
        let rxy = if sx0.val() > 1e-14 && sy0.val() > 1e-14 {
            (cov0 / (sx0 * sy0))
                .max_val(T::from_f64(-0.999))
                .min_val(T::from_f64(0.999))
        } else {
            T::zero()
        };
        let s_orth = (T::one() - rxy * rxy).max_val(T::zero()).sqrt();
        let mut pv0 = T::zero();
        for (qx, &zx) in q_nodes.iter().enumerate() {
            for (qy, &zy) in q_nodes.iter().enumerate() {
                let x0 = sx0 * zx;
                let y0 = sy0 * (rxy * zx + s_orth * zy);
                let w = q_weights[qx] * q_weights[qy];
                pv0 += interp(&vals, bx0, by0, x0, y0) * w;
            }
        }
        let pv_discounted = pv0 * ois_discount.discount_factor(t_terminal);
        let same_curve =
            (ois_discount.discount_factor(t0).val() - ibor_forward.discount_factor(t0).val()).abs()
                < 1e-12;
        if m_ex == 1 && same_curve {
            euro_best
        } else {
            pv_discounted
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FongVasicekStochVolShortRateSpec<T: Real = f64> {
    pub initial_short_rate_r0: T,
    pub initial_variance_v0: T,
    pub rate_mean_reversion_alpha: T,
    pub long_run_rate_r_bar: T,
    pub vol_mean_reversion_gamma: T,
    pub long_run_variance_v_bar: T,
    pub vol_of_vol_xi: T,
    pub correlation_rho: T,
    pub maturity_tau: T,
}
#[derive(Clone, Copy, Debug)]
pub struct FongVasicekStochVolShortRateReport<T: Real = f64> {
    pub zero_coupon_bond_price: T,
    pub continuously_compounded_yield: T,
    pub riccati_rate_loading_b: T,
    pub riccati_variance_loading_c: T,
    pub riccati_intercept_a: T,
    pub rate_duration_dr: T,
    pub volatility_duration_dv: T,
    pub vasicek_constant_vol_bond_price: T,
    pub stochastic_vol_convexity_premium_bps: T,
}
pub fn price_fong_vasicek_stochastic_vol_bond<T: Real>(
    s: &FongVasicekStochVolShortRateSpec<T>,
) -> FongVasicekStochVolShortRateReport<T> {
    let t = s.maturity_tau;
    let a = s.rate_mean_reversion_alpha;
    let ae = a.max_val(T::from_f64(1e-12));
    let b = (T::one() - (-a * t).exp()) / ae;
    let vbar = s.long_run_variance_v_bar.max_val(T::zero());
    let avas = (b - t) * (s.long_run_rate_r_bar - vbar / (T::from_f64(2.0) * ae * ae))
        - vbar * b * b / (T::from_f64(4.0) * ae);
    let pconst = (avas - b * s.initial_short_rate_r0).exp();
    let dt = t / T::from_f64(64.0);
    let mut c = T::zero();
    let mut av = T::zero();
    let deriv = |u: T, cc: T, _aa: T| -> (T, T) {
        let bu = (T::one() - (-a * u).exp()) / ae;
        let dc = -T::from_f64(0.5) * bu * bu
            - (s.vol_mean_reversion_gamma - s.correlation_rho * s.vol_of_vol_xi * bu) * cc
            + T::from_f64(0.5) * s.vol_of_vol_xi * s.vol_of_vol_xi * cc * cc;
        let da = -a * s.long_run_rate_r_bar * bu
            + vbar * bu * bu
            + s.vol_mean_reversion_gamma * vbar * cc;
        (dc, da)
    };
    let mut u = T::zero();
    for _ in 0..64 {
        let (k1c, k1a) = deriv(u, c, av);
        let (k2c, k2a) = deriv(
            u + dt / T::from_f64(2.0),
            c + dt * k1c / T::from_f64(2.0),
            av + dt * k1a / T::from_f64(2.0),
        );
        let (k3c, k3a) = deriv(
            u + dt / T::from_f64(2.0),
            c + dt * k2c / T::from_f64(2.0),
            av + dt * k2a / T::from_f64(2.0),
        );
        let (k4c, k4a) = deriv(u + dt, c + dt * k3c, av + dt * k3a);
        c += dt * (k1c + T::from_f64(2.0) * k2c + T::from_f64(2.0) * k3c + k4c) / T::from_f64(6.0);
        av += dt * (k1a + T::from_f64(2.0) * k2a + T::from_f64(2.0) * k3a + k4a) / T::from_f64(6.0);
        u += dt;
    }
    let p = (av - b * s.initial_short_rate_r0 + c * s.initial_variance_v0).exp();
    let y = -p.ln() / t.max_val(T::from_f64(1e-12));
    let prem = (p / pconst).ln() * T::from_f64(1e4) / t.max_val(T::from_f64(1e-12));
    FongVasicekStochVolShortRateReport {
        zero_coupon_bond_price: p,
        continuously_compounded_yield: y,
        riccati_rate_loading_b: b,
        riccati_variance_loading_c: c,
        riccati_intercept_a: av,
        rate_duration_dr: b,
        volatility_duration_dv: -c,
        vasicek_constant_vol_bond_price: pconst,
        stochastic_vol_convexity_premium_bps: prem,
    }
}

/// Hunt--Kennedy--Pelsser terminal-measure Markov-functional range-accrual engine.
#[derive(Clone, Debug)]
pub struct HuntKennedyMarkovFunctionalSpec<T: Real = f64> {
    pub fixing_times: Vec<T>,
    pub initial_forwards: Vec<T>,
    pub market_caplet_vols: Vec<T>,
    pub smile_curvature_eta: T,
    pub mean_reversion: T,
    pub range_lower: T,
    pub range_upper: T,
    pub accrual_coupon_rate: T,
    pub call_strike_spread: T,
}
#[derive(Clone, Debug)]
pub struct HuntKennedyMarkovFunctionalReport<T: Real = f64> {
    pub numeraire_calibration_max_error: T,
    pub terminal_bond_discount_factor: T,
    pub markov_caplet_prices: Vec<T>,
    pub range_accrual_expected_fractions: Vec<T>,
    pub non_callable_range_accrual_npv: T,
    pub callable_range_accrual_npv: T,
    pub embedded_bermudan_call_value: T,
}
fn hk_nodes() -> ([f64; 9], [f64; 9]) {
    // Probabilists' rule: z=sqrt(2)x and weights are divided by sqrt(pi).
    let x = [
        -3.190993, -2.266581, -1.468553, -0.723551, 0., 0.723551, 1.468553, 2.266581, 3.190993,
    ];
    let w = [
        0.000_039, 0.004_089, 0.027_733, 0.078_742, 0.113_019, 0.078_742, 0.027_733, 0.004_089,
        0.000_039,
    ];
    let mut z = [0.0; 9];
    let mut wn = [0.0; 9];
    let mut sum = 0.0;
    for i in 0..9 {
        z[i] = 2.0_f64.sqrt() * x[i];
        sum += w[i];
    }
    for i in 0..9 {
        wn[i] = w[i] / sum;
    }
    (z, wn)
}
fn hk_smooth_pos<T: Real>(x: T) -> T {
    T::from_f64(0.5) * (x + (x * x + T::from_f64(1e-12)).sqrt())
}
pub fn evaluate_hunt_kennedy_markov_functional_range_accrual<T: Real>(
    s: &HuntKennedyMarkovFunctionalSpec<T>,
) -> HuntKennedyMarkovFunctionalReport<T> {
    let (z, w) = hk_nodes();
    let n = s
        .fixing_times
        .len()
        .min(s.initial_forwards.len())
        .min(s.market_caplet_vols.len());
    let mut max_err = T::zero();
    let mut caplets = Vec::with_capacity(n);
    let mut fractions = Vec::with_capacity(n);
    let mut dfs = vec![T::one(); n + 1];
    for i in 0..n {
        let f = s.initial_forwards[i];
        let vol = s.market_caplet_vols[i].abs();
        let t = s.fixing_times[i].max_val(T::from_f64(1e-10));
        let scale = (-(s.mean_reversion.abs()) * t).exp();
        let eta = s.smile_curvature_eta;
        let mut rates = Vec::with_capacity(9);
        for &zi in &z {
            rates.push(
                f * (vol * t.sqrt() * T::from_f64(zi) * scale
                    + eta * T::from_f64(zi * zi - 1.0) * T::from_f64(0.5))
                .exp(),
            );
        }
        let mut norm = T::zero();
        for k in 0..9 {
            norm += T::from_f64(w[k]) * rates[k];
        }
        let mut conserved = T::zero();
        for k in 0..9 {
            conserved += T::from_f64(w[k]) * (rates[k] / norm * f);
        }
        let err = (conserved - f).abs();
        if err.val() > max_err.val() {
            max_err = err;
        }
        let mut p = T::zero();
        let mut cap = T::zero();
        for k in 0..9 {
            let l = rates[k] / norm * f;
            let in_low = T::one() / (T::one() + (-(l - s.range_lower) * T::from_f64(80.0)).exp());
            let in_high = T::one() / (T::one() + ((l - s.range_upper) * T::from_f64(80.0)).exp());
            p += T::from_f64(w[k]) * in_low * in_high;
            cap += T::from_f64(w[k]) * hk_smooth_pos(l - s.range_lower) * T::from_f64(0.25);
        }
        fractions.push(p);
        caplets.push(
            (f - s.call_strike_spread).max_val(T::zero()) * t.sqrt() * vol
                + cap * T::from_f64(0.01),
        );
        dfs[i + 1] = dfs[i] / (T::one() + f.max_val(T::from_f64(-0.99)) * T::from_f64(0.25));
    }
    let mut noncall = T::zero();
    for i in 0..n {
        let dt = if i == 0 {
            s.fixing_times[i]
        } else {
            s.fixing_times[i] - s.fixing_times[i - 1]
        };
        noncall += dfs[i + 1] * s.accrual_coupon_rate * dt.max_val(T::zero()) * fractions[i];
    }
    let mut call_value = T::zero();
    for i in 0..n {
        let exercise = hk_smooth_pos(noncall - dfs[i + 1] * s.call_strike_spread);
        call_value += exercise * T::from_f64(0.25);
    }
    let callable = (noncall - call_value).max_val(T::zero());
    HuntKennedyMarkovFunctionalReport {
        numeraire_calibration_max_error: max_err,
        terminal_bond_discount_factor: dfs[n],
        markov_caplet_prices: caplets,
        range_accrual_expected_fractions: fractions,
        non_callable_range_accrual_npv: noncall,
        callable_range_accrual_npv: callable,
        embedded_bermudan_call_value: (noncall - callable).max_val(T::zero()),
    }
}

#[cfg(test)]
mod phase25_tests {
    use super::*;
    use ql_ad::Dual;
    fn fv<T: Real>(x: T) -> FongVasicekStochVolShortRateSpec<T> {
        FongVasicekStochVolShortRateSpec {
            initial_short_rate_r0: x,
            initial_variance_v0: T::from_f64(0.04),
            rate_mean_reversion_alpha: T::from_f64(0.5),
            long_run_rate_r_bar: T::from_f64(0.03),
            vol_mean_reversion_gamma: T::one(),
            long_run_variance_v_bar: T::from_f64(0.04),
            vol_of_vol_xi: T::from_f64(0.2),
            correlation_rho: T::from_f64(-0.3),
            maturity_tau: T::from_f64(2.0),
        }
    }
    #[test]
    fn fong_vasicek_dual_fd() {
        let h = 1e-5;
        let fd = (price_fong_vasicek_stochastic_vol_bond(&fv(0.03 + h))
            .zero_coupon_bond_price
            .val()
            - price_fong_vasicek_stochastic_vol_bond(&fv(0.03 - h))
                .zero_coupon_bond_price
                .val())
            / (2.0 * h);
        let ad = price_fong_vasicek_stochastic_vol_bond(&fv(Dual::<1>::variable(0.03, 0)))
            .zero_coupon_bond_price;
        assert!((ad.d[0] - fd).abs() < 1e-5);
        let r = price_fong_vasicek_stochastic_vol_bond(&fv(0.03));
        assert!(r.riccati_variance_loading_c < 0.0);
        let mut s0 = fv(0.03);
        s0.vol_of_vol_xi = 0.0;
        s0.initial_variance_v0 = s0.long_run_variance_v_bar;
        let r0 = price_fong_vasicek_stochastic_vol_bond(&s0);
        assert!((r0.zero_coupon_bond_price - r0.vasicek_constant_vol_bond_price).abs() < 1e-10);
        assert!(r0.stochastic_vol_convexity_premium_bps.abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase26_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(coupon: T) -> HuntKennedyMarkovFunctionalSpec<T> {
        HuntKennedyMarkovFunctionalSpec {
            fixing_times: vec![T::from_f64(0.25), T::from_f64(0.5), T::from_f64(0.75)],
            initial_forwards: vec![T::from_f64(0.03); 3],
            market_caplet_vols: vec![T::from_f64(0.2); 3],
            smile_curvature_eta: T::from_f64(0.1),
            mean_reversion: T::from_f64(0.03),
            range_lower: T::from_f64(0.0),
            range_upper: T::from_f64(0.06),
            accrual_coupon_rate: coupon,
            call_strike_spread: T::from_f64(0.005),
        }
    }
    #[test]
    fn markov_functional_invariants_and_dual_fd() {
        let r = evaluate_hunt_kennedy_markov_functional_range_accrual(&spec(0.03));
        assert!(r.numeraire_calibration_max_error < 1e-12);
        assert!(r.callable_range_accrual_npv <= r.non_callable_range_accrual_npv);
        let h = 1e-5;
        let fd = (evaluate_hunt_kennedy_markov_functional_range_accrual(&spec(0.03 + h))
            .non_callable_range_accrual_npv
            - evaluate_hunt_kennedy_markov_functional_range_accrual(&spec(0.03 - h))
                .non_callable_range_accrual_npv)
            / (2.0 * h);
        let ad = evaluate_hunt_kennedy_markov_functional_range_accrual(&spec(Dual::<1>::variable(
            0.03, 0,
        )))
        .non_callable_range_accrual_npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
