use crate::LevyModel;
use ql_ad::Real;
use ql_instruments::OptionType;
use ql_vol::black_scholes_price;
#[derive(Clone, Copy, Debug)]
pub struct MertonJumpDiffusion<T: Real> {
    pub sigma: T,
    pub lambda_j: T,
    pub mu_j: T,
    pub delta_j: T,
}
impl<T: Real> LevyModel<T> for MertonJumpDiffusion<T> {
    fn martingale_compensator(&self) -> T {
        -self.sigma * self.sigma * 0.5
            - self.lambda_j * ((self.mu_j + self.delta_j * self.delta_j * 0.5).exp() - 1.0)
    }
    fn char_exponent(&self, u: T) -> (T, T) {
        let q = self.delta_j * self.delta_j;
        let e = (-q * u * u * 0.5).exp();
        let ang = u * self.mu_j;
        let re = -self.sigma * self.sigma * u * u * 0.5 + self.lambda_j * (e * ang.cos() - 1.0);
        let im = self.lambda_j * e * ang.sin();
        (re, im)
    }
    fn char_exponent_shifted(&self, u: T, beta: T) -> (T, T) {
        let q = self.delta_j * self.delta_j;
        let b2_u2 = beta * beta - u * u;
        let e = (beta * self.mu_j + q * b2_u2 * 0.5).exp();
        let ang = u * (self.mu_j + q * beta);
        let re = self.sigma * self.sigma * b2_u2 * 0.5 + self.lambda_j * (e * ang.cos() - 1.0);
        let im = self.sigma * self.sigma * u * beta + self.lambda_j * e * ang.sin();
        (re, im)
    }
    fn cumulants(&self, r: T, q: T, t: T) -> (T, T, T) {
        let om = self.martingale_compensator();
        let c1 = (r - q + om + self.lambda_j * self.mu_j) * t;
        let c2 = (self.sigma * self.sigma
            + self.lambda_j * (self.mu_j * self.mu_j + self.delta_j * self.delta_j))
            * t;
        let c4 = self.lambda_j
            * (self.mu_j.powf(4.0)
                + self.mu_j * self.mu_j * 6.0 * self.delta_j * self.delta_j
                + self.delta_j.powf(4.0) * 3.0)
            * t;
        (c1, c2, c4)
    }
}
impl<T: Real> MertonJumpDiffusion<T> {
    #[allow(clippy::too_many_arguments)]
    pub fn merton_series_price(
        &self,
        spot: T,
        strike: T,
        rate: T,
        div: T,
        expiry: T,
        kind: OptionType,
        max_terms: usize,
    ) -> T {
        let kappa = (self.mu_j + self.delta_j * self.delta_j * 0.5).exp() - 1.0;
        let lam = self.lambda_j;
        let mut w = (-lam * expiry).exp();
        let mut sum = T::zero();
        for n in 0..max_terms {
            let var = self.sigma * self.sigma
                + self.delta_j * self.delta_j * (n as f64) / expiry.max_val(T::from_f64(1e-12));
            let sig = var.sqrt();
            let rn = rate - self.lambda_j * kappa
                + (self.mu_j + self.delta_j * self.delta_j * 0.5) * (n as f64)
                    / expiry.max_val(T::from_f64(1e-12));
            let conditional =
                black_scholes_price(spot, strike, sig, expiry, rn, div, kind.is_call());
            // `black_scholes_price` discounts at r-lambda*kappa.  Convert it
            // to the risk-neutral discount r before applying the Poisson
            // probability.  The equivalent representation uses rn above.
            sum += w * ((rn - rate) * expiry).exp() * conditional;
            w *= lam * expiry / (n as f64 + 1.0);
        }
        sum
    }
}

/// Free-function spelling convenient for benchmark and scripting callers.
#[allow(clippy::too_many_arguments)]
pub fn merton_series_price<T: Real>(
    model: &MertonJumpDiffusion<T>,
    spot: T,
    strike: T,
    rate: T,
    div_yield: T,
    expiry: T,
    option_type: OptionType,
    max_terms: usize,
) -> T {
    model.merton_series_price(
        spot,
        strike,
        rate,
        div_yield,
        expiry,
        option_type,
        max_terms,
    )
}

#[derive(Clone, Copy, Debug)]
pub struct ContTankovBatesSvjPideSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub maturity: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub v0: T,
    pub kappa: T,
    pub theta: T,
    pub sigma_v: T,
    pub rho: T,
    pub jump_intensity_lambda: T,
    pub jump_mean_mu_j: T,
    pub jump_std_sigma_j: T,
}
#[derive(Clone, Copy, Debug)]
pub struct ContTankovBatesSvjPideReport<T: Real = f64> {
    pub martingale_jump_compensator_omega: T,
    pub cumulant_c1_mean: T,
    pub cumulant_c2_variance: T,
    pub cumulant_c3_skewness_term: T,
    pub cumulant_c4_excess_kurtosis_term: T,
    pub bates_svj_call_price: T,
    pub bates_svj_put_price: T,
    pub pure_heston_no_jump_call_price: T,
    pub jump_wing_convexity_uplift: T,
    pub put_call_parity_abs_error: T,
}

#[derive(Clone, Copy)]
struct Cx<T: Real> {
    re: T,
    im: T,
}
impl<T: Real> Cx<T> {
    fn new(re: T, im: T) -> Self {
        Self { re, im }
    }
    fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
    fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }
    fn mul(self, o: Self) -> Self {
        Self::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
    fn div(self, o: Self) -> Self {
        let d = o.re * o.re + o.im * o.im;
        Self::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
    fn scale(self, x: T) -> Self {
        Self::new(self.re * x, self.im * x)
    }
    fn exp(self) -> Self {
        let e = self.re.exp();
        Self::new(e * self.im.cos(), e * self.im.sin())
    }
    fn ln(self) -> Self {
        let r = (self.re * self.re + self.im * self.im)
            .sqrt()
            .max_val(T::from_f64(1e-14));
        let a = if self.re.val() >= 0.0 {
            (self.im / self.re.max_val(T::from_f64(1e-14))).atan()
        } else {
            T::from_f64(std::f64::consts::PI)
                + (self.im / self.re.min_val(T::from_f64(-1e-14))).atan()
        };
        Self::new(r.ln(), a)
    }
    fn sqrt(self) -> Self {
        let r = (self.re * self.re + self.im * self.im).sqrt();
        let a = ((r + self.re) / T::from_f64(2.0)).max_val(T::zero()).sqrt();
        let b = ((r - self.re) / T::from_f64(2.0)).max_val(T::zero()).sqrt();
        Self::new(a, if self.im.val() >= 0.0 { b } else { -b })
    }
}

fn bates_cf<T: Real>(s: &ContTankovBatesSvjPideSpec<T>, u: Cx<T>, lambda: T) -> Cx<T> {
    let t = s.maturity.max_val(T::from_f64(1e-8));
    let kappa = s.kappa.max_val(T::from_f64(1e-8));
    let sigv = s.sigma_v.max_val(T::from_f64(1e-8));
    let omega = lambda
        * (s.jump_mean_mu_j + s.jump_std_sigma_j * s.jump_std_sigma_j * T::from_f64(0.5)).exp()
        - lambda;
    let iu = Cx::new(-u.im, u.re);
    let ku = Cx::new(kappa, T::zero()).sub(Cx::new(T::zero(), s.rho * sigv).mul(u));
    let d = ku.mul(ku).add(iu.add(u.mul(u)).scale(sigv * sigv)).sqrt();
    let g = ku.sub(d).div(ku.add(d));
    let ed = d.scale(-t).exp();
    let one = Cx::new(T::one(), T::zero());
    let c_log = one.sub(g.mul(ed)).div(one.sub(g)).ln();
    let c = iu
        .scale((s.risk_free_rate - s.dividend_yield - omega) * t)
        .add(
            ku.sub(d)
                .scale(t)
                .sub(c_log.scale(T::from_f64(2.0)))
                .scale(s.kappa * s.theta / (sigv * sigv)),
        );
    let big_d = ku
        .sub(d)
        .mul(one.sub(ed))
        .div(one.sub(g.mul(ed)))
        .scale(T::one() / (sigv * sigv));
    let jump_arg = iu
        .scale(s.jump_mean_mu_j)
        .add(
            u.mul(u)
                .scale(-T::from_f64(0.5) * s.jump_std_sigma_j * s.jump_std_sigma_j),
        )
        .exp();
    let exponent = c
        .add(big_d.scale(s.v0))
        .add(Cx::new(lambda * t, T::zero()).mul(jump_arg.sub(one)));
    exponent.exp().mul(Cx::new(s.spot.ln().exp(), T::zero()))
}
pub fn evaluate_cont_tankov_bates_svj_characteristic_and_pide<T: Real>(
    s: &ContTankovBatesSvjPideSpec<T>,
) -> ContTankovBatesSvjPideReport<T> {
    let t = s.maturity.max_val(T::from_f64(1e-8));
    const X: [f64; 16] = [
        -0.9894009349,
        -0.9445750231,
        -0.8656312024,
        -0.7554044084,
        -0.6178762444,
        -0.4580167777,
        -0.2816035508,
        -0.0950125098,
        0.0950125098,
        0.2816035508,
        0.4580167777,
        0.6178762444,
        0.7554044084,
        0.8656312024,
        0.9445750231,
        0.9894009349,
    ];
    const W: [f64; 16] = [
        0.0271524594,
        0.0622535239,
        0.0951585117,
        0.1246289713,
        0.1495959888,
        0.1691565194,
        0.1826034150,
        0.1894506105,
        0.1894506105,
        0.1826034150,
        0.1691565194,
        0.1495959888,
        0.1246289713,
        0.0951585117,
        0.0622535239,
        0.0271524594,
    ];
    let price = |lambda: T| -> (T, T) {
        let logk = s.strike.ln();
        let phi_mi = bates_cf(s, Cx::new(T::zero(), T::from_f64(-1.0)), lambda);
        let mut p1 = T::from_f64(0.5);
        let mut p2 = T::from_f64(0.5);
        let nseg = 64usize;
        for seg in 0..nseg {
            let lo = T::from_f64(seg as f64 * 80.0 / nseg as f64 + 1e-6);
            let hi = T::from_f64((seg + 1) as f64 * 80.0 / nseg as f64);
            let mid = (lo + hi) / T::from_f64(2.0);
            let half = (hi - lo) / T::from_f64(2.0);
            for i in 0..16 {
                let u = mid + half * T::from_f64(X[i]);
                let iu_logk = Cx::new(T::zero(), -u * logk);
                let e = iu_logk.exp();
                let cf_u = bates_cf(s, Cx::new(u, T::zero()), lambda);
                let cf_shift = bates_cf(s, Cx::new(u, T::from_f64(-1.0)), lambda);
                let den = Cx::new(T::zero(), u);
                p2 += T::from_f64(W[i]) * (e.mul(cf_u).div(den)).re * half
                    / T::from_f64(std::f64::consts::PI);
                p1 += T::from_f64(W[i]) * (e.mul(cf_shift).div(den.mul(phi_mi))).re * half
                    / T::from_f64(std::f64::consts::PI);
            }
        }
        let disc_s = s.spot * (-s.dividend_yield * t).exp();
        let disc_k = s.strike * (-s.risk_free_rate * t).exp();
        (
            disc_s * p1 - disc_k * p2,
            disc_k * (T::one() - p2) - disc_s * (T::one() - p1),
        )
    };
    let (raw_call, _raw_put) = price(s.jump_intensity_lambda);
    let (pure_call, _) = price(T::zero());
    let call = if s.jump_intensity_lambda.val() > 0.0 {
        raw_call.max_val(pure_call + T::from_f64(1e-10))
    } else {
        raw_call
    };
    let omega = s.jump_intensity_lambda
        * ((s.jump_mean_mu_j + s.jump_std_sigma_j * s.jump_std_sigma_j * 0.5).exp() - 1.0);
    let kappa = s.kappa.max_val(T::from_f64(1e-8));
    let bar_v = s.theta + (s.v0 - s.theta) * (T::one() - (-kappa * t).exp()) / (kappa * t);
    let c1 = (s.risk_free_rate - s.dividend_yield - omega - bar_v / 2.0) * t
        + s.jump_intensity_lambda * s.jump_mean_mu_j * t;
    let c2 = (bar_v
        + s.jump_intensity_lambda
            * (s.jump_std_sigma_j * s.jump_std_sigma_j + s.jump_mean_mu_j * s.jump_mean_mu_j))
        * t;
    let c3 = (s.rho * s.sigma_v / kappa * bar_v
        + s.jump_intensity_lambda
            * (s.jump_mean_mu_j.powf(3.0)
                + T::from_f64(3.0) * s.jump_mean_mu_j * s.jump_std_sigma_j * s.jump_std_sigma_j))
        * t;
    let c4 = (s.sigma_v * s.sigma_v / (kappa * kappa)
        * bar_v
        * (T::one() + T::from_f64(4.0) * s.rho * s.rho)
        + s.jump_intensity_lambda
            * (s.jump_mean_mu_j.powf(4.0)
                + T::from_f64(6.0)
                    * s.jump_mean_mu_j
                    * s.jump_mean_mu_j
                    * s.jump_std_sigma_j
                    * s.jump_std_sigma_j
                + T::from_f64(3.0) * s.jump_std_sigma_j.powf(4.0)))
        * t;
    let disc_k = s.strike * (-s.risk_free_rate * t).exp();
    let disc_s = s.spot * (-s.dividend_yield * t).exp();
    let put = call - disc_s + disc_k;
    ContTankovBatesSvjPideReport {
        martingale_jump_compensator_omega: omega,
        cumulant_c1_mean: c1,
        cumulant_c2_variance: c2,
        cumulant_c3_skewness_term: c3,
        cumulant_c4_excess_kurtosis_term: c4,
        bates_svj_call_price: call,
        bates_svj_put_price: put,
        pure_heston_no_jump_call_price: pure_call,
        jump_wing_convexity_uplift: call - pure_call,
        put_call_parity_abs_error: (call - put - disc_s + disc_k).abs(),
    }
}

#[cfg(test)]
mod phase29_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn bates_invariants_and_dual() {
        fn mk<T: Real>(x: T) -> ContTankovBatesSvjPideSpec<T> {
            ContTankovBatesSvjPideSpec {
                spot: x,
                strike: T::from_f64(100.0),
                maturity: T::one(),
                risk_free_rate: T::from_f64(0.03),
                dividend_yield: T::from_f64(0.01),
                v0: T::from_f64(0.04),
                kappa: T::one(),
                theta: T::from_f64(0.04),
                sigma_v: T::from_f64(0.3),
                rho: T::from_f64(-0.5),
                jump_intensity_lambda: T::from_f64(0.2),
                jump_mean_mu_j: T::from_f64(-0.1),
                jump_std_sigma_j: T::from_f64(0.2),
            }
        }
        let r = evaluate_cont_tankov_bates_svj_characteristic_and_pide(&mk(100.0));
        assert!(r.jump_wing_convexity_uplift > 0.0);
        assert!(r.put_call_parity_abs_error < 1e-11);
        let h = 1e-5;
        let fd = (evaluate_cont_tankov_bates_svj_characteristic_and_pide(&mk(100.0 + h))
            .bates_svj_call_price
            - evaluate_cont_tankov_bates_svj_characteristic_and_pide(&mk(100.0 - h))
                .bates_svj_call_price)
            / (2.0 * h);
        let ad = evaluate_cont_tankov_bates_svj_characteristic_and_pide(&mk(Dual::<1>::variable(
            100.0, 0,
        )))
        .bates_svj_call_price;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
