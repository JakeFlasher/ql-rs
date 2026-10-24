use ql_ad::Real;
use ql_instruments::OptionType;

#[allow(clippy::too_many_arguments)]
pub fn margrabe_exchange_option<T: Real>(
    f1: T,
    f2: T,
    vol1: T,
    vol2: T,
    rho: T,
    rate: T,
    expiry: T,
    q1_weight: T,
    q2_weight: T,
) -> T {
    let sig = (vol1 * vol1 + vol2 * vol2 - rho * vol1 * vol2 * 2.0)
        .max_val(T::zero())
        .sqrt();
    if sig.val() <= 1e-14 || expiry.val() <= 0.0 {
        return (-rate * expiry).exp() * (q1_weight * f1 - q2_weight * f2).max_val(T::zero());
    }
    let a = q1_weight * f1;
    let b = q2_weight * f2;
    let st = sig * expiry.sqrt();
    let d1 = (a / b).ln() / st + st * 0.5;
    let d2 = d1 - st;
    (-rate * expiry).exp() * (a * d1.norm_cdf() - b * d2.norm_cdf())
}
#[allow(clippy::too_many_arguments)]
pub fn kirk_spread_option<T: Real>(
    f1: T,
    f2: T,
    strike: T,
    vol1: T,
    vol2: T,
    rho: T,
    rate: T,
    expiry: T,
    option_type: OptionType,
) -> T {
    let den = f2 + strike;
    let w = f2 / den;
    let sig = (vol1 * vol1 + w * w * vol2 * vol2 - rho * w * vol1 * vol2 * 2.0)
        .max_val(T::zero())
        .sqrt();
    if sig.val() <= 1e-14 || expiry.val() <= 0.0 {
        return (-rate * expiry).exp()
            * (if option_type.is_call() {
                f1 - den
            } else {
                den - f1
            })
            .max_val(T::zero());
    }
    let st = sig * expiry.sqrt();
    let d1 = (f1 / den).ln() / st + st * 0.5;
    let d2 = d1 - st;
    let s = if option_type.is_call() { 1.0 } else { -1.0 };
    (-rate * expiry).exp()
        * T::from_f64(s)
        * (f1 * (T::from_f64(s) * d1).norm_cdf() - den * (T::from_f64(s) * d2).norm_cdf())
}
#[allow(clippy::too_many_arguments)]
pub fn spark_spread_option<T: Real>(
    power: T,
    gas: T,
    heat_rate: T,
    strike: T,
    vol_power: T,
    vol_gas: T,
    rho: T,
    rate: T,
    expiry: T,
    option_type: OptionType,
) -> T {
    kirk_spread_option(
        power,
        gas * heat_rate,
        strike,
        vol_power,
        vol_gas,
        rho,
        rate,
        expiry,
        option_type,
    )
}

#[derive(Clone, Debug)]
pub struct ThreeLegSpreadOptionSpec<T: Real> {
    pub forward1: T,
    pub forward2: T,
    pub forward3: T,
    pub weight1: T,
    pub weight2: T,
    pub weight3: T,
    pub strike: T,
    pub vol1: T,
    pub vol2: T,
    pub vol3: T,
    pub rho12: T,
    pub rho13: T,
    pub rho23: T,
    pub rate: T,
    pub expiry: T,
    pub option_type: OptionType,
}
pub fn price_three_leg_spread_option<T: Real>(s: &ThreeLegSpreadOptionSpec<T>) -> T {
    if s.weight3.val().abs() <= 1e-14 {
        return kirk_spread_option(
            s.forward1 * s.weight1,
            s.forward2 * s.weight2,
            s.strike,
            s.vol1,
            s.vol2,
            s.rho12,
            s.rate,
            s.expiry,
            s.option_type,
        );
    }
    let a = s.forward1 * s.weight1 + s.forward3 * s.weight3;
    let b = s.forward2 * s.weight2 + s.strike;
    let va = (s.forward1 * s.weight1 / a).powf(2.0) * s.vol1 * s.vol1
        + (s.forward3 * s.weight3 / a).powf(2.0) * s.vol3 * s.vol3
        + T::from_f64(2.0) * s.forward1 * s.weight1 * s.forward3 * s.weight3 / (a * a)
            * s.rho13
            * s.vol1
            * s.vol3;
    let vb = s.forward2 * s.weight2 / b * s.forward2 * s.weight2 / b * s.vol2 * s.vol2;
    let cov = (s.forward1 * s.weight1 / a)
        * (s.forward2 * s.weight2 / b)
        * s.rho12
        * s.vol1
        * s.vol2
        + (s.forward3 * s.weight3 / a) * (s.forward2 * s.weight2 / b) * s.rho23 * s.vol3 * s.vol2;
    let sig = (va + vb - T::from_f64(2.0) * cov).max_val(T::zero()).sqrt();
    let st = sig * s.expiry.sqrt();
    if st.val() <= 1e-14 {
        return (-s.rate * s.expiry).exp()
            * (if s.option_type.is_call() {
                a - b
            } else {
                b - a
            })
            .max_val(T::zero());
    }
    let d1 = (a / b).ln() / st + st * 0.5;
    let d2 = d1 - st;
    let sign = T::from_f64(if s.option_type.is_call() { 1.0 } else { -1.0 });
    (-s.rate * s.expiry).exp() * sign * (a * (sign * d1).norm_cdf() - b * (sign * d2).norm_cdf())
}
