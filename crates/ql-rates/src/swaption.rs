use ql_ad::Real;
use ql_curves::YieldTermStructure;

fn normal_sabr_vol<T: Real>(f: T, k: T, expiry: f64, alpha: T, beta: f64, rho: T, nu: T) -> T {
    let fk = f * k;
    let fb = fk.powf((1.0 - beta) * 0.5);
    let z = nu / alpha * fb * (f - k);
    let ratio = if z.val().abs() < 1e-7 {
        T::one() - rho * z * 0.5 + (T::from_f64(2.0) - rho * rho * 3.0) * z * z / 12.0
    } else {
        let disc = (T::one() - rho * z * 2.0 + z * z).sqrt();
        z / ((disc + z - rho) / (T::one() - rho)).ln()
    };
    let omb = T::from_f64(1.0 - beta);
    let corr = T::one()
        + T::from_f64(expiry)
            * (omb * omb * alpha * alpha / (fk.pow(omb) * 24.0)
                + rho * T::from_f64(beta) * nu * alpha / (fk.pow(omb * 0.5) * 4.0)
                + (T::from_f64(2.0) - rho * rho * 3.0) * nu * nu / 24.0);
    alpha * fk.powf(beta * 0.5) * ratio * corr
}
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SwaptionSettlement {
    Physical,
    CashParYield,
}
#[derive(Clone, Debug)]
pub struct EuropeanSwaption {
    pub notional: f64,
    pub strike: f64,
    pub expiry: f64,
    pub swap_times: Vec<f64>,
    pub is_payer: bool,
    pub settlement: SwaptionSettlement,
}
impl EuropeanSwaption {
    fn forward_annuity<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        ois: &D,
        fwd: &F,
    ) -> (T, T, T) {
        if self.swap_times.len() < 2 {
            return (T::zero(), T::zero(), T::zero());
        }
        let mut a = T::zero();
        let mut num = T::zero();
        for w in self.swap_times.windows(2) {
            let tau = w[1] - w[0];
            let ff = (fwd.discount_factor(w[0]) / fwd.discount_factor(w[1]) - T::one()) / tau;
            let df = ois.discount_factor(w[1]);
            a += T::from_f64(tau) * df;
            num += T::from_f64(tau) * ff * df;
        }
        let s = num / a;
        let ann = match self.settlement {
            SwaptionSettlement::Physical => a,
            SwaptionSettlement::CashParYield => {
                let mut x = T::zero();
                for (i, w) in self.swap_times.windows(2).enumerate() {
                    let tau = w[1] - w[0];
                    x += ois.discount_factor(self.swap_times[0]) * T::from_f64(tau)
                        / (T::one() + T::from_f64(tau) * s).powf((i + 1) as f64);
                }
                x
            }
        };
        (s, ann, num)
    }
    pub fn price_bachelier<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        ois: &D,
        fwd: &F,
        normal_vol: T,
    ) -> T {
        let (s, a, _) = self.forward_annuity(ois, fwd);
        let st = normal_vol * self.expiry.sqrt();
        let om = if self.is_payer { 1.0 } else { -1.0 };
        if st.val().abs() <= 1e-14 {
            return T::from_f64(self.notional)
                * a
                * (T::from_f64(om) * (s - self.strike)).max_val(T::zero());
        }
        let d = if st.val().abs() > 1e-14 {
            (s - self.strike) / st
        } else {
            T::zero()
        };
        T::from_f64(self.notional)
            * a
            * (T::from_f64(om) * (s - self.strike) * (T::from_f64(om) * d).norm_cdf()
                + st * d.norm_pdf())
    }
    pub fn price_shifted_black76<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        ois: &D,
        fwd: &F,
        vol: T,
        shift: f64,
    ) -> T {
        let (s, a, _) = self.forward_annuity(ois, fwd);
        let st = vol * self.expiry.sqrt();
        let f = s + shift;
        let k = T::from_f64(self.strike + shift);
        if st.val().abs() <= 1e-14 || f.val() <= 0.0 || k.val() <= 0.0 {
            let om = if self.is_payer { 1.0 } else { -1.0 };
            return T::from_f64(self.notional)
                * a
                * (T::from_f64(om) * (s - self.strike)).max_val(T::zero());
        }
        let d1 = ((f / k).ln() + T::from_f64(0.5) * st * st) / st;
        let d2 = d1 - st;
        let om = if self.is_payer { 1.0 } else { -1.0 };
        T::from_f64(self.notional)
            * a
            * (T::from_f64(om)
                * (f * (T::from_f64(om) * d1).norm_cdf() - k * (T::from_f64(om) * d2).norm_cdf()))
    }
    pub fn price_sabr_bachelier<T: Real, D: YieldTermStructure<T>, F: YieldTermStructure<T>>(
        &self,
        ois: &D,
        fwd: &F,
        alpha: T,
        beta: f64,
        rho: T,
        nu: T,
    ) -> T {
        let (s, _, _) = self.forward_annuity(ois, fwd);
        let vol = normal_sabr_vol(
            s,
            T::from_f64(self.strike),
            self.expiry,
            alpha,
            beta,
            rho,
            nu,
        );
        self.price_bachelier(ois, fwd, vol)
    }
}
