#![allow(
    clippy::needless_range_loop,
    clippy::assign_op_pattern,
    clippy::approx_constant
)]
use ql_ad::Real;

#[derive(Clone, Copy, Debug)]
pub struct GatheralJacquierSviDensitySpec<T: Real = f64> {
    pub forward_price: T,
    pub strike_price: T,
    pub expiry_years: T,
    pub svi_a: T,
    pub svi_b: T,
    pub svi_rho: T,
    pub svi_m: T,
    pub svi_sigma: T,
}

#[derive(Clone, Copy, Debug)]
pub struct GatheralJacquierSviDensityReport<T: Real = f64> {
    pub log_moneyness_k: T,
    pub svi_total_variance_w: T,
    pub svi_implied_vol: T,
    pub durrleman_butterfly_g: T,
    pub breeden_litzenberger_pdf: T,
    pub roger_lee_right_wing_margin: T,
    pub roger_lee_left_wing_margin: T,
}

pub fn evaluate_gatheral_jacquier_arbitrage_free_svi_density<T: Real>(
    spec: &GatheralJacquierSviDensitySpec<T>,
) -> GatheralJacquierSviDensityReport<T> {
    let f = spec.forward_price.max_val(T::from_f64(1e-6));
    let k_str = spec.strike_price.max_val(T::from_f64(1e-6));
    let mat = spec.expiry_years.max_val(T::from_f64(1e-6));
    let k = (k_str / f).ln();
    let km = k - spec.svi_m;
    let sig = spec.svi_sigma.max_val(T::from_f64(1e-4));
    let rad = (km * km + sig * sig).sqrt();

    let w = (spec.svi_a + spec.svi_b * (spec.svi_rho * km + rad)).max_val(T::from_f64(1e-6));
    let wp = spec.svi_b * (spec.svi_rho + km / rad);
    let wpp = spec.svi_b * sig * sig / (rad * rad * rad);

    let term1_root = T::one() - (k * wp) / (T::from_f64(2.0) * w);
    let term1 = term1_root * term1_root;
    let term2 = T::from_f64(0.25) * wp * wp * (T::one() / w + T::from_f64(0.25));
    let g = term1 - term2 + T::from_f64(0.5) * wpp;

    let sqrt_w = w.sqrt();
    let d2 = -k / sqrt_w - T::from_f64(0.5) * sqrt_w;
    let inv_sqrt_2pi = T::from_f64(0.398_942_280_401_432_7);
    let pdf = (g.max_val(T::from_f64(1e-8)) / sqrt_w)
        * inv_sqrt_2pi
        * (-T::from_f64(0.5) * d2 * d2).exp();

    let abs_rho = (spec.svi_rho * spec.svi_rho + T::from_f64(1e-16)).sqrt();
    let right_margin = T::from_f64(2.0) - spec.svi_b * (T::one() + abs_rho);
    let left_margin = T::from_f64(2.0) - spec.svi_b * (T::one() - spec.svi_rho);

    GatheralJacquierSviDensityReport {
        log_moneyness_k: k,
        svi_total_variance_w: w,
        svi_implied_vol: (w / mat).sqrt(),
        durrleman_butterfly_g: g,
        breeden_litzenberger_pdf: pdf,
        roger_lee_right_wing_margin: right_margin,
        roger_lee_left_wing_margin: left_margin,
    }
}

#[cfg(test)]
mod phase38_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gatheral_jacquier_svi_density_invariants_and_ad_parity() {
        let spec = GatheralJacquierSviDensitySpec {
            forward_price: 100.0,
            strike_price: 102.0,
            expiry_years: 1.0,
            svi_a: 0.015,
            svi_b: 0.12,
            svi_rho: -0.35,
            svi_m: 0.02,
            svi_sigma: 0.22,
        };
        let rep = evaluate_gatheral_jacquier_arbitrage_free_svi_density(&spec);
        assert!(rep.durrleman_butterfly_g > 0.0);
        assert!(rep.breeden_litzenberger_pdf > 0.0);
        assert!(rep.roger_lee_right_wing_margin > 0.0);

        let h = 1e-6;
        let r_ad = evaluate_gatheral_jacquier_arbitrage_free_svi_density(
            &GatheralJacquierSviDensitySpec {
                forward_price: Dual::<1>::from_value(spec.forward_price),
                strike_price: Dual::<1>::variable(spec.strike_price, 0),
                expiry_years: Dual::<1>::from_value(spec.expiry_years),
                svi_a: Dual::<1>::from_value(spec.svi_a),
                svi_b: Dual::<1>::from_value(spec.svi_b),
                svi_rho: Dual::<1>::from_value(spec.svi_rho),
                svi_m: Dual::<1>::from_value(spec.svi_m),
                svi_sigma: Dual::<1>::from_value(spec.svi_sigma),
            },
        );
        let mut sp = spec;
        sp.strike_price += h;
        let mut sm = spec;
        sm.strike_price -= h;
        let fd = (evaluate_gatheral_jacquier_arbitrage_free_svi_density(&sp)
            .breeden_litzenberger_pdf
            - evaluate_gatheral_jacquier_arbitrage_free_svi_density(&sm).breeden_litzenberger_pdf)
            / (2.0 * h);
        assert!((r_ad.breeden_litzenberger_pdf.d[0] - fd).abs() < 1e-5);
    }
}
