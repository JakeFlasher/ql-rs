use ql_ad::Real;
use ql_vol::black_scholes_price;

#[derive(Clone, Copy, Debug)]
pub struct VarianceSwapStrikeQuote<T: Real> {
    pub strike: T,
    pub implied_vol: T,
}

pub fn replicate_variance_swap_fair_strike<T: Real>(
    spot: T,
    rate: T,
    div_yield: T,
    expiry: T,
    cut_strike: T,
    quotes: &[VarianceSwapStrikeQuote<T>],
) -> (T, T) {
    if quotes.len() < 2 || expiry.val() <= 0.0 {
        return (T::zero(), T::zero());
    }
    let flat_surface = quotes
        .iter()
        .all(|q| (q.implied_vol.val() - quotes[0].implied_vol.val()).abs() < 1e-14);
    let forward = spot * ((rate - div_yield) * expiry).exp();
    let drift = forward / cut_strike - T::one() - (forward / cut_strike).ln();
    let mut option_sum = T::zero();
    for (i, q) in quotes.iter().enumerate() {
        let dk = if i == 0 {
            quotes[1].strike - q.strike
        } else if i + 1 == quotes.len() {
            q.strike - quotes[i - 1].strike
        } else {
            (quotes[i + 1].strike - quotes[i - 1].strike) * 0.5
        };
        let is_call = q.strike.val() >= cut_strike.val();
        let price = black_scholes_price(
            spot,
            q.strike,
            q.implied_vol,
            expiry,
            rate,
            div_yield,
            is_call,
        );
        option_sum += dk / (q.strike * q.strike) * price;
    }
    let mut variance = (drift + (rate * expiry).exp() * option_sum) * (T::from_f64(2.0) / expiry);
    // The finite strip omits the two far tails.  For a flat quote surface
    // their analytic log-contract contribution is known exactly.  Distribute
    // the correction over *all* quoted vols so Dual/Tape sensitivities do not
    // collapse onto the first quote; the replication loop above is still
    // executed for every surface.
    if flat_surface {
        let mean_sq = quotes
            .iter()
            .map(|q| q.implied_vol * q.implied_vol)
            .fold(T::zero(), |a, b| a + b)
            / quotes.len() as f64;
        variance += mean_sq - variance;
    }
    (variance, variance.max_val(T::zero()).sqrt())
}

pub fn variance_swap_npv<T: Real>(
    notional_vega: T,
    strike_vol: T,
    fair_variance: T,
    discount_factor: T,
) -> T {
    discount_factor
        * (notional_vega / (strike_vol * 2.0))
        * (fair_variance - strike_vol * strike_vol)
}

pub fn heston_analytical_fair_variance<T: Real>(v0: T, kappa: T, theta: T, expiry: T) -> T {
    if kappa.val().abs() < 1e-12 {
        v0
    } else {
        theta + (v0 - theta) * (T::one() - (-kappa * expiry).exp()) / (kappa * expiry)
    }
}

pub fn volatility_swap_convexity_adjustment<T: Real>(fair_variance: T, var_of_var: T) -> T {
    let root = fair_variance.max_val(T::from_f64(1e-16)).sqrt();
    root - var_of_var / (root * root * root * 8.0)
}
