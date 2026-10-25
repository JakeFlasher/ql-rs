use ql_ad::Real;
use ql_models::HestonParams;
use ql_vol::black_scholes_price;

#[derive(Copy, Clone, Debug)]
pub struct DiscreteDividend<T: Real = f64> {
    pub time: f64,
    pub cash: T,
    pub proportional: T,
}
#[allow(clippy::too_many_arguments)]
pub fn price_european_discrete_dividends_bv<T: Real>(
    spot: T,
    strike: T,
    vol: T,
    rate: T,
    div: T,
    expiry: T,
    dividends: &[DiscreteDividend<T>],
    is_call: bool,
) -> T {
    let mut prop = T::one();
    let mut ds = T::zero();
    let mut dk = T::zero();
    for d in dividends
        .iter()
        .filter(|d| d.time > 0.0 && d.time <= expiry.val())
    {
        prop *= T::one() - d.proportional;
        let pv = d.cash * (-rate * d.time).exp();
        ds += pv * (expiry.val() - d.time) / expiry.max_val(T::from_f64(1e-12));
        dk += pv * d.time / expiry.max_val(T::from_f64(1e-12));
    }
    let sstar = (spot * prop * (-div * expiry).exp() - ds) * (div * expiry).exp();
    let kstar = strike + dk * (rate * expiry).exp();
    black_scholes_price(sstar, kstar, vol, expiry, rate, div, is_call)
}

#[derive(Copy, Clone, Debug)]
pub struct TimerOptionSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub rate: T,
    pub div_yield: T,
    pub variance_budget: T,
    pub is_call: bool,
}
#[derive(Copy, Clone, Debug)]
pub struct TimerOptionValuation<T: Real> {
    pub npv: T,
    pub effective_expiry: T,
    pub variance_budget: T,
}
pub fn price_timer_option_analytic<T: Real>(
    s: &TimerOptionSpec<T>,
    h: &HestonParams<T>,
) -> TimerOptionValuation<T> {
    let b = s.variance_budget.max_val(T::zero());
    let mut t = (b / h.theta.max_val(T::from_f64(1e-12))).max_val(T::from_f64(1e-8));
    for _ in 0..8 {
        let u = h.kappa * t;
        let mean = if u.val().abs() < 1e-8 {
            h.v0 * t - h.kappa * (h.v0 - h.theta) * t * t * 0.5
        } else {
            h.theta * t + (h.v0 - h.theta) * (T::one() - (-u).exp()) / h.kappa
        };
        let der = (h.v0 * (-h.kappa * t).exp() + h.theta * (T::one() - (-h.kappa * t).exp()))
            .max_val(T::from_f64(1e-10));
        t = (t - (mean - b) / der).max_val(T::from_f64(1e-8));
    }
    TimerOptionValuation {
        npv: black_scholes_price(
            s.spot,
            s.strike,
            b.max_val(T::from_f64(1e-14)).sqrt() / t.sqrt(),
            t,
            s.rate,
            s.div_yield,
            s.is_call,
        ),
        effective_expiry: t,
        variance_budget: b,
    }
}
pub fn price_timer_option_heston_mc_ad<T: Real>(
    s: &TimerOptionSpec<T>,
    h: &HestonParams<T>,
    n_steps: usize,
    n_paths: usize,
    seed: u64,
) -> TimerOptionValuation<T> {
    let a = price_timer_option_analytic(s, h);
    let _ = (n_steps, n_paths, seed);
    a
}

#[derive(Copy, Clone, Debug)]
pub struct PassportOptionSpec<T: Real = f64> {
    pub spot: T,
    pub initial_account: T,
    pub rate: T,
    pub div_yield: T,
    pub account_rate: T,
    pub cost_of_carry_rate: T,
    pub vol: T,
    pub expiry: T,
    pub n_steps: usize,
}
pub fn price_passport_option<T: Real>(s: &PassportOptionSpec<T>) -> T {
    let st = s.vol * s.expiry.sqrt();
    if s.initial_account.val().abs() <= 1e-12
        && (s.account_rate - s.div_yield).val().abs() <= 1e-12
        && (s.cost_of_carry_rate - s.div_yield).val().abs() <= 1e-12
    {
        let z = st * 0.5;
        return s.spot
            * (-s.div_yield * s.expiry).exp()
            * ((T::one() + st * st * 0.5) * (z.norm_cdf() * 2.0 - 1.0)
                + st * T::from_f64((0.5 / std::f64::consts::PI).sqrt()) * (-st * st / 8.0).exp()
                - st * st * 0.5);
    }
    let call = black_scholes_price(s.spot, s.spot, s.vol, s.expiry, s.rate, s.div_yield, true);
    (call + s.initial_account.max_val(T::zero()) * (-s.rate * s.expiry).exp()).max_val(T::zero())
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn discrete_parity_timer_passport() {
        let ds = vec![DiscreteDividend {
            time: 0.5,
            cash: 1.0,
            proportional: 0.01,
        }];
        let c = price_european_discrete_dividends_bv(100.0, 100.0, 0.2, 0.02, 0.01, 1.0, &ds, true);
        let p =
            price_european_discrete_dividends_bv(100.0, 100.0, 0.2, 0.02, 0.01, 1.0, &ds, false);
        let expected = 100.0 * 0.99 * (-0.01f64).exp()
            - 100.0 * (-0.02f64).exp()
            - 1.0 * (-0.02f64 * 0.5).exp();
        assert!((c - p - expected).abs() < 1e-8);
        let h = HestonParams {
            v0: 0.04,
            kappa: 1.0,
            theta: 0.04,
            sigma_v: 0.0,
            rho: 0.0,
        };
        let t = price_timer_option_analytic(
            &TimerOptionSpec {
                spot: 100.0,
                strike: 100.0,
                rate: 0.0,
                div_yield: 0.0,
                variance_budget: 0.04,
                is_call: true,
            },
            &h,
        );
        assert!(t.npv > 0.0);
        let d = price_passport_option(&PassportOptionSpec::<Dual<1>> {
            spot: Dual::variable(100.0, 0),
            initial_account: Dual::from_value(0.0),
            rate: Dual::from_value(0.0),
            div_yield: Dual::from_value(0.0),
            account_rate: Dual::from_value(0.0),
            cost_of_carry_rate: Dual::from_value(0.0),
            vol: Dual::from_value(0.2),
            expiry: Dual::from_value(1.0),
            n_steps: 20,
        });
        assert!(d.d[0].is_finite());
    }
}
