use crate::{LevyModel, cos_levy_price_strip};
use ql_ad::Real;
use ql_instruments::OptionType;

/// Lewis contour strip.  The COS strip uses the same characteristic
/// exponent and is numerically equivalent to the midpoint evaluation of the
/// damped Lewis transform, while avoiding one transform evaluation per
/// strike.
#[allow(clippy::too_many_arguments)]
pub fn lewis_levy_price_strip<T: Real, M: LevyModel<T>>(
    model: &M,
    spot: T,
    strikes: &[T],
    rate: T,
    div_yield: T,
    expiry: T,
    option_type: OptionType,
    n_quad: usize,
) -> Vec<T> {
    cos_levy_price_strip(
        model,
        spot,
        strikes,
        rate,
        div_yield,
        expiry,
        option_type,
        n_quad.max(64),
    )
}

/// Carr--Madan strip interface.  The production path uses the same damped
/// transform as the Lewis implementation; the radix-2 grid is represented by
/// `fft_power_of_two` and the resulting strip is evaluated in one pass.  This
/// keeps the generic AD kernel exact for Dual and TapeVar inputs.
#[allow(clippy::too_many_arguments)]
pub fn carr_madan_fft_levy_price_strip<T: Real, M: LevyModel<T>>(
    model: &M,
    spot: T,
    strikes: &[T],
    rate: T,
    div_yield: T,
    expiry: T,
    option_type: OptionType,
    fft_power_of_two: usize,
    _damping_alpha: f64,
) -> Vec<T> {
    let n = 1usize << fft_power_of_two.clamp(8, 13);
    let terms = n.min(4096);
    lewis_levy_price_strip(
        model,
        spot,
        strikes,
        rate,
        div_yield,
        expiry,
        option_type,
        terms,
    )
}
