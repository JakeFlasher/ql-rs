use ql_ad::{Real, ift_lift_1d};
use ql_math::brent_solve;
use ql_mc::PhiloxRng;
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
fn erlang_cdf<T: Real>(x: T, mu: T, n: usize) -> T {
    if n == 0 {
        return T::one();
    }
    let y = (mu * x).max_val(T::zero());
    T::one()
        - (-y).exp()
            * (0..n)
                .map(|k| y.powf(k as f64) / fact(k))
                .fold(T::zero(), |a, b| a + b)
}
fn erlang_tail<T: Real>(x: T, mu: T, n: usize) -> T {
    T::one() - erlang_cdf(x, mu, n)
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
        self.payment_times
            .iter()
            .enumerate()
            .map(|(i, &t)| {
                let a = T::from_f64(*self.accrual_fractions.get(i).unwrap_or(&1.0));
                let red = if i + 1 == self.payment_times.len() {
                    T::one()
                } else {
                    T::zero()
                };
                self.face_amount * (self.coupon_rate * a + red) * (-self.discount_rate * t).exp()
            })
            .fold(T::zero(), |a, b| a + b)
    }
    pub fn price_analytical(&self) -> CatBondValuation<T> {
        if matches!(self.peril, CatPerilSpec::PoissonPareto { .. }) {
            return self.price_mc(4096, 0x9e3779b97f4a7c15);
        }
        let (lambda, ms) = match self.peril {
            CatPerilSpec::PoissonExponential {
                lambda,
                mean_severity,
            } => (lambda, mean_severity),
            _ => unreachable!(),
        };
        let mu = T::one() / ms.max_val(T::from_f64(1e-12));
        let mut price = T::zero();
        let mut loss = T::zero();
        let mut attach = T::zero();
        let mut exhaust = T::zero();
        for (i, &t) in self.payment_times.iter().enumerate() {
            let a = T::from_f64(*self.accrual_fractions.get(i).unwrap_or(&1.0));
            let mut surv = T::zero();
            let mut ap = T::zero();
            let mut ep = T::zero();
            for n in 0..=256 {
                let w = poisson_weight(lambda, t, n);
                let (at, ex) = match self.risk {
                    NotionalRisk::Digital { threshold } => (threshold, threshold),
                    NotionalRisk::Proportional {
                        attachment,
                        exhaustion,
                    } => (attachment, exhaustion),
                };
                let pa = erlang_tail(at, mu, n);
                let pe = erlang_tail(ex, mu, n);
                ap += w * pa;
                ep += w * pe;
                let rem = match self.risk {
                    NotionalRisk::Digital { threshold } => erlang_cdf(threshold, mu, n),
                    NotionalRisk::Proportional { .. } => {
                        if n == 0 {
                            T::one()
                        } else {
                            let ea =
                                (T::from_f64(n as f64) / mu) * erlang_tail(at, mu, n + 1) - at * pa;
                            let ee =
                                (T::from_f64(n as f64) / mu) * erlang_tail(ex, mu, n + 1) - ex * pe;
                            T::one() - (ea - ee) / (ex - at).max_val(T::from_f64(1e-12))
                        }
                    }
                };
                surv += w * rem;
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
    pub fn price_mc(&self, n_paths: usize, seed: u64) -> CatBondValuation<T> {
        let np = n_paths.max(1);
        let rng = PhiloxRng::new(seed);
        let mut p = 0.0;
        let mut l = 0.0;
        let mut a_prob = 0.0;
        let mut e_prob = 0.0;
        for path in 0..np {
            let mut sev = 0.0;
            let mut te = 0.0;
            let mut step = 0u32;
            for (i, &t) in self.payment_times.iter().enumerate() {
                loop {
                    let (lam, sev_kind) = match self.peril {
                        CatPerilSpec::PoissonExponential { lambda, .. } => (lambda.val(), 0),
                        CatPerilSpec::PoissonPareto { lambda, .. } => (lambda.val(), 1),
                    };
                    if lam <= 0.0 {
                        break;
                    }
                    let u = rng.uniform4(path as u64, step)[0].clamp(1e-15, 1. - 1e-15);
                    step += 1;
                    te += -u.ln() / lam;
                    if te > t {
                        break;
                    }
                    let u2 = rng.uniform4(path as u64, step)[1].clamp(1e-15, 1. - 1e-15);
                    step += 1;
                    sev += if sev_kind == 0 {
                        match self.peril {
                            CatPerilSpec::PoissonExponential { mean_severity, .. } => {
                                -mean_severity.val() * u2.ln()
                            }
                            _ => 0.,
                        }
                    } else {
                        match self.peril {
                            CatPerilSpec::PoissonPareto {
                                scale_xm,
                                shape_alpha,
                                ..
                            } => scale_xm.val() * u2.powf(-1. / shape_alpha.val().max(1e-12)),
                            _ => 0.,
                        }
                    };
                }
                let (surv, at, ex) = match self.risk {
                    NotionalRisk::Digital { threshold } => (
                        if sev <= threshold.val() { 1. } else { 0. },
                        sev > threshold.val(),
                        sev > threshold.val(),
                    ),
                    NotionalRisk::Proportional {
                        attachment,
                        exhaustion,
                    } => {
                        let aa = attachment.val();
                        let ee = exhaustion.val();
                        (
                            (1. - (sev - aa).max(0.) / (ee - aa).max(1e-12)).clamp(0., 1.),
                            sev > aa,
                            sev > ee,
                        )
                    }
                };
                let df = (-self.discount_rate.val() * t).exp();
                let accr = *self.accrual_fractions.get(i).unwrap_or(&1.);
                let red = if i + 1 == self.payment_times.len() {
                    surv
                } else {
                    0.
                };
                p += self.face_amount.val() * (self.coupon_rate.val() * accr * surv + red) * df;
                l += self.face_amount.val() * (1. - surv) * df;
                a_prob += at as i32 as f64;
                e_prob += ex as i32 as f64;
            }
        }
        let den = np as f64;
        let denp = den * self.payment_times.len().max(1) as f64;
        CatBondValuation {
            price: T::from_f64(p / den),
            attachment_prob: T::from_f64(a_prob / denp),
            exhaustion_prob: T::from_f64(e_prob / denp),
            expected_loss: T::from_f64(l / den),
        }
    }
}
pub fn price_analytical<T: Real>(s: &CatBondSpec<T>) -> CatBondValuation<T> {
    s.price_analytical()
}
pub fn price_mc<T: Real>(s: &CatBondSpec<T>, n: usize, seed: u64) -> CatBondValuation<T> {
    s.price_mc(n, seed)
}
pub fn cat_bond_z_spread<T: Real>(s: &CatBondSpec<T>, market: T) -> Result<T, String> {
    let root = brent_solve(
        |z| {
            let mut sx = s.clone();
            sx.discount_rate = s.discount_rate + T::from_f64(z);
            sx.price_analytical().price.val() - market.val()
        },
        -1.,
        1.,
        1e-12,
        200,
    )
    .map_err(|e| e.to_string())?;
    let mut sx = s.clone();
    sx.discount_rate = s.discount_rate + T::from_f64(root);
    let mut fx = T::zero();
    let mut fxx = 0.;
    for (i, &t) in s.payment_times.iter().enumerate() {
        let a = T::from_f64(*s.accrual_fractions.get(i).unwrap_or(&1.));
        let red = if i + 1 == s.payment_times.len() {
            T::one()
        } else {
            T::zero()
        };
        let cf = s.face_amount * (s.coupon_rate * a + red);
        let df = (-sx.discount_rate * T::from_f64(t)).exp();
        fx -= T::from_f64(t) * cf * df;
        fxx += t * t * cf.val() * df.val();
    }
    Ok(ift_lift_1d(
        root,
        |z| {
            let mut q = s.clone();
            q.discount_rate = s.discount_rate + z;
            q.price_analytical().price - market
        },
        fx,
        fxx,
    ))
}
