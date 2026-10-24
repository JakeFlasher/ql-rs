use crate::BatesParams;
use ql_ad::Real;
pub fn replicate_cboe_vix_index<T: Real>(expiry: T, rate: T, forward: T, quotes: &[(T, T)]) -> T {
    // CBOE's strip replication requires a positive expiry, forward and strike
    // grid.  Returning a zero quote for an unusable strip is preferable to
    // propagating NaNs through a larger calibration graph.
    if quotes.len() < 2
        || !expiry.val().is_finite()
        || expiry.val() <= 0.0
        || !forward.val().is_finite()
        || forward.val() <= 0.0
        || !rate.val().is_finite()
    {
        return T::zero();
    }
    let mut q = quotes.to_vec();
    q.sort_by(|a, b| a.0.val().total_cmp(&b.0.val()));
    if q.iter()
        .any(|(k, p)| !k.val().is_finite() || k.val() <= 0.0 || !p.val().is_finite())
    {
        return T::zero();
    }
    // Duplicate strikes do not define a finite integration cell.  Keep the
    // first occurrence to avoid double-counting their option quote.
    q.dedup_by(|a, b| (a.0.val() - b.0.val()).abs() <= 1e-14);
    if q.len() < 2 {
        return T::zero();
    }
    let mut sum = T::zero();
    for i in 0..q.len() {
        let k = q[i].0;
        let dk = if i == 0 {
            q[1].0 - k
        } else if i + 1 == q.len() {
            k - q[i - 1].0
        } else {
            (q[i + 1].0 - q[i - 1].0) * 0.5
        };
        if dk.val() <= 0.0 {
            return T::zero();
        }
        sum += dk / (k * k) * q[i].1.max_val(T::zero());
    }
    let k0 = q
        .iter()
        .filter(|(k, _)| k.val() <= forward.val())
        .map(|(k, _)| *k)
        .last()
        .unwrap_or(q[0].0);
    let var = (sum * T::from_f64(2.0) / expiry) * (rate * expiry).exp()
        - (forward / k0 - T::one()).powf(2.0) / expiry;
    T::from_f64(100.0) * var.max_val(T::zero()).sqrt()
}
#[derive(Copy, Clone, Debug)]
pub struct VixOptionValuation<T: Real> {
    pub vix_future: T,
    pub call_price: T,
    pub put_price: T,
}
pub fn heston_bates_vix_future_and_option_price<T: Real>(
    p: &BatesParams<T>,
    expiry: T,
    vix_window: f64,
    strike: T,
    discount_factor: T,
) -> VixOptionValuation<T> {
    if !expiry.val().is_finite()
        || expiry.val() < 0.0
        || !vix_window.is_finite()
        || vix_window <= 0.0
        || !strike.val().is_finite()
        || strike.val() <= 0.0
        || !discount_factor.val().is_finite()
        || discount_factor.val() < 0.0
    {
        return VixOptionValuation {
            vix_future: T::zero(),
            call_price: T::zero(),
            put_price: T::zero(),
        };
    }
    let k = p.heston.kappa;
    let a = if k.val().abs() < 1e-8 {
        T::one()
    } else {
        (T::one() - (-k * T::from_f64(vix_window)).exp()) / (k * T::from_f64(vix_window))
    };
    let jump = p.lambda * (p.mu_j.exp() - T::one() - p.mu_j + p.sigma_j * p.sigma_j * 0.5);
    let b = p.heston.theta * (T::one() - a) + T::from_f64(2.0) * jump;
    let mean =
        p.heston.v0 * (-k * expiry).exp() + p.heston.theta * (T::one() - (-k * expiry).exp());
    let f = (a * mean + b).max_val(T::zero()).sqrt();
    if expiry.val() <= 0.0 {
        return VixOptionValuation {
            vix_future: f,
            call_price: (f - strike).max_val(T::zero()) * discount_factor,
            put_price: (strike - f).max_val(T::zero()) * discount_factor,
        };
    }
    let vvol = (p.heston.sigma_v.abs() * expiry.sqrt() * 0.5).max_val(T::from_f64(1e-6));
    let c = ql_vol::black_scholes_price(f, strike, vvol, expiry, T::zero(), T::zero(), true)
        * discount_factor;
    let put = ql_vol::black_scholes_price(f, strike, vvol, expiry, T::zero(), T::zero(), false)
        * discount_factor;
    VixOptionValuation {
        vix_future: f,
        call_price: c,
        put_price: put,
    }
}
