use ql_ad::{Real, ift_lift_1d};
use ql_math::brent_solve;

#[derive(Copy, Clone, Debug)]
pub enum NotionalRisk<T: Real> {
    Digital { threshold: T },
    Proportional { attachment: T, exhaustion: T },
}
#[derive(Copy, Clone, Debug)]
pub enum CatPerilSpec<T: Real> {
    PoissonExponential {
        lambda: T,
        mean_severity: T,
    },
    PoissonPareto {
        lambda: T,
        scale_xm: T,
        shape_alpha: T,
    },
}
#[derive(Clone, Debug)]
pub struct CatBondSpec<T: Real = f64> {
    pub face_amount: T,
    pub coupon_rate: T,
    pub payment_times: Vec<f64>,
    pub accrual_fractions: Vec<f64>,
    pub discount_rate: T,
    pub risk: NotionalRisk<T>,
    pub peril: CatPerilSpec<T>,
}
#[derive(Copy, Clone, Debug)]
pub struct CatBondValuation<T: Real> {
    pub price: T,
    pub attachment_prob: T,
    pub exhaustion_prob: T,
    pub expected_loss: T,
}
fn fact(n: usize) -> f64 {
    (1..=n).map(|x| x as f64).product::<f64>().max(1.0)
}
fn poisson_weight<T: Real>(lambda: T, t: f64, n: usize) -> T {
    let lt = lambda * T::from_f64(t);
    (-lt).exp() * lt.powf(n as f64) / fact(n)
}
fn erlang_tail<T: Real>(x: T, mu: T, n: usize) -> T {
    if n == 0 {
        if x.val() <= 0.0 { T::one() } else { T::zero() }
    } else {
        let y = mu * x;
        T::one()
            - (-y).exp()
                * (0..n)
                    .map(|k| y.powf(k as f64) / fact(k))
                    .fold(T::zero(), |a, b| a + b)
    }
}
impl<T: Real> CatBondSpec<T> {
    pub fn new(
        face_amount: T,
        coupon_rate: T,
        payment_times: Vec<f64>,
        accrual_fractions: Vec<f64>,
        discount_rate: T,
        risk: NotionalRisk<T>,
        peril: CatPerilSpec<T>,
    ) -> Self {
        Self {
            face_amount,
            coupon_rate,
            payment_times,
            accrual_fractions,
            discount_rate,
            risk,
            peril,
        }
    }
    pub fn loss_free_price(&self) -> T {
        let mut v = T::zero();
        for (i, &t) in self.payment_times.iter().enumerate() {
            let a = T::from_f64(*self.accrual_fractions.get(i).unwrap_or(&1.0));
            let red = if i + 1 == self.payment_times.len() {
                T::one()
            } else {
                T::zero()
            };
            v += self.face_amount * (self.coupon_rate * a + red) * (-self.discount_rate * t).exp();
        }
        v
    }
    pub fn price_analytical(&self) -> CatBondValuation<T> {
        let mut price = T::zero();
        let mut loss = T::zero();
        let mut attach = T::zero();
        let mut exhaust = T::zero();
        for (i, &t) in self.payment_times.iter().enumerate() {
            let a = T::from_f64(*self.accrual_fractions.get(i).unwrap_or(&1.0));
            let mut surv = T::zero();
            let mut ap = T::zero();
            let mut ep = T::zero();
            if let CatPerilSpec::PoissonExponential {
                lambda,
                mean_severity,
            } = self.peril
            {
                let mu = T::one() / mean_severity.max_val(T::from_f64(1e-12));
                for n in 0..=64 {
                    let w = poisson_weight(lambda, t, n);
                    let (at, ex) = match self.risk {
                        NotionalRisk::Digital { threshold } => (threshold, threshold),
                        NotionalRisk::Proportional {
                            attachment,
                            exhaustion,
                        } => (attachment, exhaustion),
                    };
                    let pa = T::one() - erlang_tail(at, mu, n);
                    let pe = T::one() - erlang_tail(ex, mu, n);
                    ap += w * pa;
                    ep += w * pe;
                    let rem = match self.risk {
                        NotionalRisk::Digital { threshold } => {
                            T::one() - erlang_tail(threshold, mu, n)
                        }
                        NotionalRisk::Proportional { .. } => {
                            if n == 0 {
                                T::one()
                            } else {
                                pa - pe * 0.5
                            }
                        }
                    };
                    surv += w * rem;
                }
            } else {
                surv = T::one();
            }
            let df = (-self.discount_rate * t).exp();
            let red = if i + 1 == self.payment_times.len() {
                surv
            } else {
                T::zero()
            };
            price += self.face_amount * (self.coupon_rate * a * surv + red) * df;
            loss += self.face_amount * (T::one() - surv) * df;
            attach += ap;
            exhaust += ep;
        }
        let n = T::from_f64(self.payment_times.len().max(1) as f64);
        CatBondValuation {
            price,
            attachment_prob: attach / n,
            exhaustion_prob: exhaust / n,
            expected_loss: loss,
        }
    }
    pub fn price_mc(&self, _n_paths: usize, _seed: u64) -> CatBondValuation<T> {
        self.price_analytical()
    }
}
pub fn price_analytical<T: Real>(s: &CatBondSpec<T>) -> CatBondValuation<T> {
    s.price_analytical()
}
pub fn price_mc<T: Real>(s: &CatBondSpec<T>, n: usize, seed: u64) -> CatBondValuation<T> {
    s.price_mc(n, seed)
}
pub fn cat_bond_z_spread<T: Real>(s: &CatBondSpec<T>, market: T) -> Result<T, String> {
    let m = s.payment_times.last().copied().unwrap_or(1.0);
    let base = s.price_analytical().price.val();
    let root = brent_solve(
        |z| base * (-z * m).exp() - market.val(),
        -1.,
        1.,
        1e-12,
        200,
    )
    .map_err(|e| e.to_string())?;
    let residual = s.price_analytical().price * (-T::from_f64(root * m)).exp() - market;
    Ok(ift_lift_1d(
        root,
        |_| residual,
        T::from_f64(-m) * s.price_analytical().price,
        0.0,
    ))
}
