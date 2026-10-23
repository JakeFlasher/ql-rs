use crate::black::VolError;
use implied_vol::{DefaultSpecialFn, ImpliedNormalVolatility};
use ql_ad::{Real, ift_lift_1d};
pub fn bachelier_price<T: Real>(
    forward: T,
    strike: T,
    normal_vol: T,
    expiry: T,
    discount: T,
    is_call: bool,
) -> T {
    let s = normal_vol * expiry.sqrt();
    if s.val() <= 0.0 {
        return discount
            * if is_call {
                (forward - strike).max_val(T::zero())
            } else {
                (strike - forward).max_val(T::zero())
            };
    }
    let d = (forward - strike) / s;
    let p = if is_call {
        (forward - strike) * d.norm_cdf() + s * d.norm_pdf()
    } else {
        (strike - forward) * (-d).norm_cdf() + s * d.norm_pdf()
    };
    discount * p
}
pub fn implied_normal_vol<T: Real>(
    forward: T,
    strike: T,
    expiry: T,
    undiscounted_price: T,
    is_call: bool,
) -> Result<T, VolError> {
    let sigma = ImpliedNormalVolatility::builder()
        .forward(forward.val())
        .strike(strike.val())
        .expiry(expiry.val())
        .option_price(undiscounted_price.val())
        .is_call(is_call)
        .build()
        .and_then(|x| x.calculate::<DefaultSpecialFn>())
        .ok_or(VolError::SolverFailed)?;
    let residual =
        |v: T| bachelier_price(forward, strike, v, expiry, T::one(), is_call) - undiscounted_price;
    let ss = T::from_f64(sigma);
    let sq = expiry.sqrt();
    let d = (forward - strike) / (ss * sq);
    let vega = sq * d.norm_pdf();
    let volga = vega * d * d / ss;
    Ok(ift_lift_1d(sigma, residual, vega, volga.val()))
}
pub fn bachelier_greeks<T: Real>(
    forward: T,
    strike: T,
    vol: T,
    expiry: T,
    discount: T,
    is_call: bool,
) -> (T, T, T, T, T, T, T, T) {
    let price = bachelier_price(forward, strike, vol, expiry, discount, is_call);
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        let sign = if is_call { 1.0 } else { -1.0 };
        let indicator = if forward.val() > strike.val() {
            1.0
        } else if forward.val() < strike.val() {
            0.0
        } else {
            0.5
        };
        return (
            price,
            discount * T::from_f64(sign * indicator),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::zero(),
        );
    }
    let s = vol * expiry.sqrt();
    let d = (forward - strike) / s;
    let sign = T::from_f64(if is_call { 1.0 } else { -1.0 });
    let delta = discount * sign * d.norm_cdf();
    let gamma = discount * d.norm_pdf() / s;
    let vega = discount * expiry.sqrt() * d.norm_pdf();
    let theta = discount * (-T::from_f64(0.5) * vol / expiry.sqrt()) * d.norm_pdf();
    let rho = T::zero();
    let vanna = T::zero();
    let volga = discount * expiry.sqrt() * d.norm_pdf() * d * d / vol;
    (price, delta, gamma, vega, theta, rho, vanna, volga)
}
