use ql_ad::Real;
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
    fn survival_exponential(&self, t: f64, lambda: T, mean_severity: T) -> (T, T, T) {
        let mu = T::one() / mean_severity.max_val(T::from_f64(1e-12));
        let (at, ex) = match self.risk {
            NotionalRisk::Digital { threshold } => (threshold, threshold),
            NotionalRisk::Proportional {
                attachment,
                exhaustion,
            } => (attachment, exhaustion),
        };
        let mut surv = T::zero();
        let mut ap = T::zero();
        let mut ep = T::zero();
        for n in 0..=256 {
            let w = poisson_weight(lambda, t, n);
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
        (surv, ap, ep)
    }
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
            } => (lambda.max_val(T::zero()), mean_severity),
            _ => unreachable!(),
        };
        let mut price = T::zero();
        let mut loss = T::zero();
        let mut attach = T::zero();
        let mut exhaust = T::zero();
        let mut prev_surv = T::one();
        for (i, &t) in self.payment_times.iter().enumerate() {
            let a = T::from_f64(*self.accrual_fractions.get(i).unwrap_or(&1.0));
            let (surv, ap, ep) = self.survival_exponential(t, lambda, ms);
            let df = (-self.discount_rate * t).exp();
            let red = if i + 1 == self.payment_times.len() {
                surv
            } else {
                T::zero()
            };
            price += self.face_amount * (self.coupon_rate * a * surv + red) * df;
            loss += self.face_amount * (prev_surv - surv) * df;
            prev_surv = surv;
            if i + 1 == self.payment_times.len() {
                attach = ap;
                exhaust = ep;
            }
        }
        CatBondValuation {
            price,
            attachment_prob: attach,
            exhaustion_prob: exhaust,
            expected_loss: loss,
        }
    }
    pub fn price_mc(&self, n_paths: usize, seed: u64) -> CatBondValuation<T> {
        let np = n_paths.max(1);
        let rng = PhiloxRng::new(seed);
        let mut p = T::zero();
        let mut l = T::zero();
        let mut a_prob = T::zero();
        let mut e_prob = T::zero();
        for path in 0..np {
            let mut sev = T::zero();
            let mut prev_surv = T::one();
            let mut step = 0u32;
            for (i, &t) in self.payment_times.iter().enumerate() {
                let prev_t = if i == 0 {
                    0.0
                } else {
                    self.payment_times[i - 1]
                };
                let mut te = T::from_f64(prev_t);
                loop {
                    let (lambda, sev_kind) = match self.peril {
                        CatPerilSpec::PoissonExponential { lambda, .. } => (lambda, 0),
                        CatPerilSpec::PoissonPareto { lambda, .. } => (lambda, 1),
                    };
                    if lambda.val() <= 0.0 {
                        break;
                    }
                    let u = rng.uniform4(path as u64, step)[0].clamp(1e-15, 1. - 1e-15);
                    step += 1;
                    te += T::from_f64(-u.ln()) / lambda;
                    if te.val() > t {
                        break;
                    }
                    let u2 = rng.uniform4(path as u64, step)[1].clamp(1e-15, 1. - 1e-15);
                    step += 1;
                    sev += if sev_kind == 0 {
                        match self.peril {
                            CatPerilSpec::PoissonExponential { mean_severity, .. } => {
                                mean_severity * T::from_f64(-u2.ln())
                            }
                            _ => T::zero(),
                        }
                    } else {
                        match self.peril {
                            CatPerilSpec::PoissonPareto {
                                scale_xm,
                                shape_alpha,
                                ..
                            } => {
                                scale_xm
                                    * T::from_f64(u2)
                                        .pow(-T::one() / shape_alpha.max_val(T::from_f64(1e-12)))
                            }
                            _ => T::zero(),
                        }
                    };
                }
                let (surv, at, ex) = match self.risk {
                    NotionalRisk::Digital { threshold } => {
                        let hit = sev.val() > threshold.val();
                        (if hit { T::zero() } else { T::one() }, hit, hit)
                    }
                    NotionalRisk::Proportional {
                        attachment,
                        exhaustion,
                    } => {
                        let aa = attachment;
                        let ee = exhaustion;
                        (
                            (T::one()
                                - (sev - aa).max_val(T::zero())
                                    / (ee - aa).max_val(T::from_f64(1e-12)))
                            .max_val(T::zero())
                            .min_val(T::one()),
                            sev.val() > aa.val(),
                            sev.val() > ee.val(),
                        )
                    }
                };
                let df = (-self.discount_rate * T::from_f64(t)).exp();
                let accr = T::from_f64(*self.accrual_fractions.get(i).unwrap_or(&1.));
                let red = if i + 1 == self.payment_times.len() {
                    surv
                } else {
                    T::zero()
                };
                p += self.face_amount * (self.coupon_rate * accr * surv + red) * df;
                l += self.face_amount * (prev_surv - surv) * df;
                prev_surv = surv;
                if i + 1 == self.payment_times.len() {
                    a_prob += T::from_f64(at as u8 as f64);
                    e_prob += T::from_f64(ex as u8 as f64);
                }
            }
        }
        let den = T::from_f64(np as f64);
        CatBondValuation {
            price: p / den,
            attachment_prob: a_prob / den,
            exhaustion_prob: e_prob / den,
            expected_loss: l / den,
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
    ql_math::brent_solve_ad(
        |z| {
            let mut sx = s.clone();
            sx.discount_rate = s.discount_rate + z;
            sx.price_analytical().price - market
        },
        -1.0,
        1.0,
        1e-12,
        200,
    )
    .map_err(|e| e.to_string())
}
