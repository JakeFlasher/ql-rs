use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_math::MathError;
#[derive(Copy, Clone, Debug)]
pub struct HullWhite1F<T: Real> {
    pub a: T,
    pub sigma: T,
}
impl<T: Real> HullWhite1F<T> {
    pub fn b_factor(&self, t: f64, big_t: f64) -> T {
        let dt = big_t - t;
        if self.a.val().abs() < 1e-12 {
            T::from_f64(dt)
        } else {
            (T::one() - (-self.a * dt).exp()) / self.a
        }
    }
    pub fn zcb_price(&self, curve: &impl YieldTermStructure<T>, t: f64, big_t: f64, r_t: T) -> T {
        let b = self.b_factor(t, big_t);
        let p0t = curve.discount_factor(t);
        let p0_t = curve.discount_factor(big_t);
        let f = curve.instantaneous_forward(t);
        let var_r = if self.a.val().abs() < 1e-10 {
            T::from_f64(t)
        } else {
            (T::one() - (-self.a * T::from_f64(2.0 * t)).exp()) / (self.a * T::from_f64(2.0))
        };
        let a = (p0_t / p0t)
            * (f * b - T::from_f64(0.5) * self.sigma * self.sigma * b * b * var_r).exp();
        a * (-b * r_t).exp()
    }
    pub fn zcb_call_put(
        &self,
        curve: &impl YieldTermStructure<T>,
        expiry: f64,
        maturity: f64,
        strike: T,
        is_call: bool,
    ) -> T {
        let pe = curve.discount_factor(expiry);
        let pm = curve.discount_factor(maturity);
        let b = self.b_factor(expiry, maturity);
        let var_r = if self.a.val().abs() < 1e-10 {
            T::from_f64(expiry)
        } else {
            (T::one() - (-self.a * T::from_f64(2.0 * expiry)).exp()) / (self.a * T::from_f64(2.0))
        };
        let sigma_p = self.sigma * b * var_r.sqrt();
        if sigma_p.val().abs() < 1e-14 {
            let intrinsic = pm - strike * pe;
            return if is_call {
                intrinsic.max_val(T::zero())
            } else {
                (-intrinsic).max_val(T::zero())
            };
        }
        let h = ((pm / (pe * strike)).ln() + T::from_f64(0.5) * sigma_p * sigma_p) / sigma_p;
        if is_call {
            pm * h.norm_cdf() - strike * pe * (h - sigma_p).norm_cdf()
        } else {
            strike * pe * (-h + sigma_p).norm_cdf() - pm * (-h).norm_cdf()
        }
    }
    pub fn swaption_jamshidian(
        &self,
        curve: &impl YieldTermStructure<T>,
        expiry: f64,
        payment_times: &[f64],
        accrual: &[f64],
        fixed: T,
        payer: bool,
    ) -> Result<T, MathError> {
        if payment_times.is_empty() {
            return Ok(T::zero());
        }
        let f = |r: T| {
            let mut sum = T::zero();
            for (i, &tm) in payment_times.iter().enumerate() {
                let tau = accrual.get(i).copied().unwrap_or(1.0);
                let ci = if i + 1 == payment_times.len() {
                    fixed * tau + T::one()
                } else {
                    fixed * tau
                };
                sum += ci * self.zcb_price(curve, expiry, tm, r);
            }
            sum - T::one()
        };
        let rstar = ql_math::brent_solve_ad(f, -0.5, 1.0, 1e-12, 200)?;
        let mut total = T::zero();
        for (i, &tm) in payment_times.iter().enumerate() {
            let tau = accrual.get(i).copied().unwrap_or(1.0);
            let ci = if i + 1 == payment_times.len() {
                fixed * tau + T::one()
            } else {
                fixed * tau
            };
            let ki = self.zcb_price(curve, expiry, tm, rstar);
            total += ci * self.zcb_call_put(curve, expiry, tm, ki, !payer);
        }
        Ok(total)
    }
}
