use implied_vol::{DefaultSpecialFn, ImpliedBlackVolatility};
use ql_ad::{Dual2, Real, ift_lift_1d};
use thiserror::Error;
#[derive(Debug, Error, Clone)]
pub enum VolError {
    #[error("invalid option inputs")]
    InvalidInput,
    #[error("implied volatility solver failed")]
    SolverFailed,
}
#[inline]
pub fn black76_undiscounted<T: Real>(forward: T, strike: T, vol: T, expiry: T, is_call: bool) -> T {
    let intrinsic = if is_call {
        (forward - strike).max_val(T::zero())
    } else {
        (strike - forward).max_val(T::zero())
    };
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        return intrinsic;
    }
    let s = expiry.sqrt();
    let vs = vol * s;
    let d1 = ((forward / strike).ln() + T::from_f64(0.5) * vol * vol * expiry) / vs;
    let d2 = d1 - vs;
    if is_call {
        forward * d1.norm_cdf() - strike * d2.norm_cdf()
    } else {
        strike * (-d2).norm_cdf() - forward * (-d1).norm_cdf()
    }
}
#[inline]
pub fn black76_price<T: Real>(
    forward: T,
    strike: T,
    vol: T,
    expiry: T,
    discount: T,
    is_call: bool,
) -> T {
    discount * black76_undiscounted(forward, strike, vol, expiry, is_call)
}
pub fn black_scholes_price<T: Real>(
    spot: T,
    strike: T,
    vol: T,
    expiry: T,
    rate: T,
    dividend: T,
    is_call: bool,
) -> T {
    let f = spot * ((rate - dividend) * expiry).exp();
    let d = (-(rate) * expiry).exp();
    black76_price(f, strike, vol, expiry, d, is_call)
}
pub fn implied_black_vol<T: Real>(
    forward: T,
    strike: T,
    expiry: T,
    undiscounted_price: T,
    is_call: bool,
) -> Result<T, VolError> {
    let sigma = ImpliedBlackVolatility::builder()
        .forward(forward.val())
        .strike(strike.val())
        .expiry(expiry.val())
        .option_price(undiscounted_price.val())
        .is_call(is_call)
        .build()
        .and_then(|x| x.calculate::<DefaultSpecialFn>())
        .ok_or(VolError::SolverFailed)?;
    let residual =
        |v: T| black76_undiscounted(forward, strike, v, expiry, is_call) - undiscounted_price;
    let ss = T::from_f64(sigma);
    let sq = expiry.sqrt();
    let d1 = ((forward / strike).ln() + T::from_f64(0.5) * ss * ss * expiry) / (ss * sq);
    let vega = forward * d1.norm_pdf() * sq;
    let volga = vega * d1 * (d1 - ss * sq) / ss;
    Ok(ift_lift_1d(sigma, residual, vega, volga.val()))
}

/// Explicit four-input second-order wrapper for benchmark and risk callers.
pub fn implied_black_vol_dual2(
    forward: Dual2<4>,
    strike: Dual2<4>,
    expiry: Dual2<4>,
    undiscounted_price: Dual2<4>,
    is_call: bool,
) -> Result<Dual2<4>, VolError> {
    implied_black_vol(forward, strike, expiry, undiscounted_price, is_call)
}

pub fn black76_greeks<T: Real>(
    forward: T,
    strike: T,
    vol: T,
    expiry: T,
    discount: T,
    is_call: bool,
) -> (T, T, T, T, T, T, T, T) {
    let p = black76_price(forward, strike, vol, expiry, discount, is_call);
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
            p,
            discount * T::from_f64(sign * indicator),
            T::zero(),
            T::zero(),
            T::zero(),
            -expiry * p,
            T::zero(),
            T::zero(),
        );
    }
    let s = expiry.sqrt();
    let d1 = ((forward / strike).ln() + T::from_f64(0.5) * vol * vol * expiry) / (vol * s);
    let d2 = d1 - vol * s;
    let pdf = d1.norm_pdf();
    let c = T::from_f64(if is_call { 1.0 } else { -1.0 });
    let delta = discount * c * d1.norm_cdf();
    let gamma = discount * pdf / (forward * vol * s);
    let vega = discount * forward * pdf * s;
    let theta = -(discount * forward * pdf * vol) / (T::from_f64(2.0) * s);
    let rho = -expiry * p;
    let vanna = discount * pdf * s * (T::one() - d1 / (vol * s));
    let volga = vega * d1 * d2 / vol;
    (p, delta, gamma, vega, theta, rho, vanna, volga)
}
