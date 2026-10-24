use crate::BatesParams;
use ql_ad::Real;
pub fn replicate_cboe_vix_index<T: Real>(expiry: T, rate: T, forward: T, quotes: &[(T, T)]) -> T {
    if quotes.len() < 2 {
        return T::zero();
    }
    let mut sum = T::zero();
    for i in 0..quotes.len() {
        let k = quotes[i].0;
        let dk = if i == 0 {
            quotes[1].0 - k
        } else if i + 1 == quotes.len() {
            k - quotes[i - 1].0
        } else {
            (quotes[i + 1].0 - quotes[i - 1].0) * 0.5
        };
        sum += dk / (k * k) * quotes[i].1;
    }
    let k0 = quotes
        .iter()
        .filter(|(k, _)| k.val() <= forward.val())
        .map(|(k, _)| *k)
        .last()
        .unwrap_or(quotes[0].0);
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
    let c = (f - strike).max_val(T::zero()) * discount_factor;
    let put = (strike - f).max_val(T::zero()) * discount_factor;
    VixOptionValuation {
        vix_future: f,
        call_price: c,
        put_price: put,
    }
}
