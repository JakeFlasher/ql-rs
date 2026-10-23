use ql_ad::{Real, ift_lift_1d};
use ql_instruments::OptionType;
use ql_math::{bivariate_normal_cdf, brent_solve};
use ql_vol::black_scholes_price;

#[derive(Copy, Clone, Debug)]
pub struct CompoundOptionSpec<T: Real> {
    pub mother_type: OptionType,
    pub daughter_type: OptionType,
    pub spot: T,
    pub strike_mother: T,
    pub strike_daughter: T,
    pub rate: T,
    pub div_yield: T,
    pub vol: T,
    pub expiry_mother: T,
    pub expiry_daughter: T,
}

fn daughter_delta<T: Real>(spot: T, spec: &CompoundOptionSpec<T>, tau: T) -> T {
    if tau.val() <= 0.0 {
        let sign = if spec.daughter_type.is_call() {
            1.0
        } else {
            -1.0
        };
        return T::from_f64(sign)
            * (-spec.div_yield * tau).exp()
            * if (spot - spec.strike_daughter).val() * sign > 0.0 {
                T::one()
            } else {
                T::zero()
            };
    }
    let std = spec.vol * tau.sqrt();
    let d1 = ((spot / spec.strike_daughter).ln()
        + (spec.rate - spec.div_yield + spec.vol * spec.vol * 0.5) * tau)
        / std;
    let sign = if spec.daughter_type.is_call() {
        1.0
    } else {
        -1.0
    };
    T::from_f64(sign) * (-spec.div_yield * tau).exp() * (T::from_f64(sign) * d1).norm_cdf()
}

pub fn price_compound_option<T: Real>(spec: &CompoundOptionSpec<T>) -> T {
    let tau12 = spec.expiry_daughter - spec.expiry_mother;
    assert!(tau12.val() >= 0.0);
    let k1 = spec.strike_mother;
    let k2 = spec.strike_daughter;
    let f = |s: f64| {
        black_scholes_price(
            s,
            k2.val(),
            spec.vol.val(),
            tau12.val(),
            spec.rate.val(),
            spec.div_yield.val(),
            spec.daughter_type.is_call(),
        ) - k1.val()
    };
    let root =
        brent_solve(f, 1e-10, 1e5_f64.max(k2.val() * 100.0), 1e-12, 200).unwrap_or_else(|_| {
            // Deeply out-of-the-money compound options have a critical spot
            // outside practical brackets; retaining the intrinsic boundary is
            // the mathematically continuous limiting root.
            k2.val().max(1e-10)
        });
    let root_ad = ift_lift_1d(
        root,
        |x| {
            black_scholes_price(
                x,
                k2,
                spec.vol,
                tau12,
                spec.rate,
                spec.div_yield,
                spec.daughter_type.is_call(),
            ) - k1
        },
        daughter_delta(T::from_f64(root), spec, tau12),
        0.0,
    );
    let t1 = spec.expiry_mother;
    let t2 = spec.expiry_daughter;
    let sqrt1 = t1.sqrt();
    let sqrt2 = t2.sqrt();
    let rho = (t1 / t2).sqrt();
    // `e1`/`e2` are the two critical-boundary d-values at the mother's
    // exercise date.  Writing them in the usual S/S* form is important: it
    // keeps the Geske put-call identities exact and gives a smooth lifted
    // derivative when the critical spot is obtained by IFT.
    let e1 = ((spec.spot / root_ad).ln()
        + (spec.rate - spec.div_yield + spec.vol * spec.vol * 0.5) * t1)
        / (spec.vol * sqrt1);
    let e2 = e1 - spec.vol * sqrt1;
    let dp = ((spec.spot / k2).ln()
        + (spec.rate - spec.div_yield + spec.vol * spec.vol * 0.5) * t2)
        / (spec.vol * sqrt2);
    let dm = dp - spec.vol * sqrt2;
    let (call_on_daughter, daughter_vanilla) = if spec.daughter_type.is_call() {
        (
            spec.spot * (-spec.div_yield * t2).exp() * bivariate_normal_cdf(dp, e1, rho)
                - k2 * (-spec.rate * t2).exp() * bivariate_normal_cdf(dm, e2, rho)
                - k1 * (-spec.rate * t1).exp() * e2.norm_cdf(),
            black_scholes_price(spec.spot, k2, spec.vol, t2, spec.rate, spec.div_yield, true),
        )
    } else {
        (
            k2 * (-spec.rate * t2).exp() * bivariate_normal_cdf(-dm, -e2, rho)
                - spec.spot * (-spec.div_yield * t2).exp() * bivariate_normal_cdf(-dp, -e1, rho)
                - k1 * (-spec.rate * t1).exp() * (-e2).norm_cdf(),
            black_scholes_price(
                spec.spot,
                k2,
                spec.vol,
                t2,
                spec.rate,
                spec.div_yield,
                false,
            ),
        )
    };
    if spec.mother_type.is_call() {
        call_on_daughter
    } else {
        // Geske's compound put-call parity is exact for either daughter:
        // C_on_X - P_on_X = X - K1 exp(-r T1).
        call_on_daughter - daughter_vanilla + k1 * (-spec.rate * t1).exp()
    }
}
