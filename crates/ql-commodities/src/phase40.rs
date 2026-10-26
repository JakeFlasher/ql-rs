//! Eydeland & Wolyniec (2003) & Geman (2005) Power/Gas Clean Spark & Tolling Spread Option
//! via Bjerksund–Stensland (2014) / Margrabe Three-Leg Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct EydelandCleanSparkSpreadOptionSpec<T: Real = f64> {
    pub forward_power_price_fp: T,
    pub forward_gas_price_fg: T,
    pub plant_heat_rate_h: T,
    pub fixed_vom_and_carbon_k: T,
    pub power_volatility_sigma_p: T,
    pub gas_volatility_sigma_g: T,
    pub power_gas_corr_rho: T,
    pub risk_free_rate: T,
    pub expiry_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct EydelandCleanSparkSpreadOptionReport<T: Real = f64> {
    pub intrinsic_clean_spark_spread: T,
    pub composite_fuel_and_vom_strike_a: T,
    pub composite_spread_volatility: T,
    pub margrabe_zero_strike_exchange_call: T,
    pub bjerksund_stensland_spark_call: T,
    pub margrabe_limit_consistency_error: T,
}

#[allow(clippy::too_many_arguments)]
fn price_spark<T: Real>(fp: T, fuel_leg: T, k_vom: T, sp: T, sg: T, rho: T, df: T, t: T) -> (T, T) {
    let a = (fuel_leg + k_vom).max_val(T::from_f64(1e-6));
    let w = fuel_leg / a;
    let var_c = (sp * sp - T::from_f64(2.0) * rho * sp * sg * w + sg * sg * w * w)
        .max_val(T::from_f64(1e-8));
    let sig_c = var_c.sqrt();
    let std = sig_c * t.sqrt();
    let d1 = ((fp / a).ln() + T::from_f64(0.5) * std * std) / std;
    let d2 = d1 - std;
    let call = df * (fp * norm_cdf(d1) - a * norm_cdf(d2));
    (call, sig_c)
}

pub fn evaluate_eydeland_clean_spark_spread_option<T: Real>(
    spec: &EydelandCleanSparkSpreadOptionSpec<T>,
) -> EydelandCleanSparkSpreadOptionReport<T> {
    let fp = spec.forward_power_price_fp.max_val(T::from_f64(1e-4));
    let fuel = (spec.plant_heat_rate_h * spec.forward_gas_price_fg).max_val(T::from_f64(1e-4));
    let t = spec.expiry_years_t.max_val(T::from_f64(0.02));
    let df = (-spec.risk_free_rate * t).exp();

    let intrinsic = fp - fuel - spec.fixed_vom_and_carbon_k;
    let a = fuel + spec.fixed_vom_and_carbon_k;

    let (bs_call, sig_comp) = price_spark(
        fp,
        fuel,
        spec.fixed_vom_and_carbon_k,
        spec.power_volatility_sigma_p,
        spec.gas_volatility_sigma_g,
        spec.power_gas_corr_rho,
        df,
        t,
    );
    let (margrabe_call, _) = price_spark(
        fp,
        fuel,
        T::zero(),
        spec.power_volatility_sigma_p,
        spec.gas_volatility_sigma_g,
        spec.power_gas_corr_rho,
        df,
        t,
    );

    // Exact Margrabe (1978) exchange formula at K_VOM = 0:
    let sig_m = (spec.power_volatility_sigma_p * spec.power_volatility_sigma_p
        - T::from_f64(2.0)
            * spec.power_gas_corr_rho
            * spec.power_volatility_sigma_p
            * spec.gas_volatility_sigma_g
        + spec.gas_volatility_sigma_g * spec.gas_volatility_sigma_g)
        .max_val(T::from_f64(1e-8))
        .sqrt();
    let std_m = sig_m * t.sqrt();
    let d1_m = ((fp / fuel).ln() + T::from_f64(0.5) * std_m * std_m) / std_m;
    let d2_m = d1_m - std_m;
    let exact_margrabe = df * (fp * norm_cdf(d1_m) - fuel * norm_cdf(d2_m));
    let m_err = (margrabe_call - exact_margrabe).abs();

    EydelandCleanSparkSpreadOptionReport {
        intrinsic_clean_spark_spread: intrinsic,
        composite_fuel_and_vom_strike_a: a,
        composite_spread_volatility: sig_comp,
        margrabe_zero_strike_exchange_call: margrabe_call,
        bjerksund_stensland_spark_call: bs_call,
        margrabe_limit_consistency_error: m_err,
    }
}

#[cfg(test)]
mod phase40_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn eydeland_clean_spark_spread_option_invariants_and_ad_parity() {
        let spec = EydelandCleanSparkSpreadOptionSpec {
            forward_power_price_fp: 68.0,
            forward_gas_price_fg: 7.2,
            plant_heat_rate_h: 7.5,
            fixed_vom_and_carbon_k: 8.5,
            power_volatility_sigma_p: 0.42,
            gas_volatility_sigma_g: 0.36,
            power_gas_corr_rho: 0.72,
            risk_free_rate: 0.04,
            expiry_years_t: 0.5,
        };
        let rep = evaluate_eydeland_clean_spark_spread_option(&spec);
        assert!(rep.margrabe_limit_consistency_error < 1e-12);
        assert!(rep.margrabe_zero_strike_exchange_call > rep.bjerksund_stensland_spark_call);

        let h = 1e-6;
        let r_ad =
            evaluate_eydeland_clean_spark_spread_option(&EydelandCleanSparkSpreadOptionSpec {
                forward_power_price_fp: Dual::<1>::variable(spec.forward_power_price_fp, 0),
                forward_gas_price_fg: Dual::<1>::from_value(spec.forward_gas_price_fg),
                plant_heat_rate_h: Dual::<1>::from_value(spec.plant_heat_rate_h),
                fixed_vom_and_carbon_k: Dual::<1>::from_value(spec.fixed_vom_and_carbon_k),
                power_volatility_sigma_p: Dual::<1>::from_value(spec.power_volatility_sigma_p),
                gas_volatility_sigma_g: Dual::<1>::from_value(spec.gas_volatility_sigma_g),
                power_gas_corr_rho: Dual::<1>::from_value(spec.power_gas_corr_rho),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                expiry_years_t: Dual::<1>::from_value(spec.expiry_years_t),
            });
        let mut sp = spec;
        sp.forward_power_price_fp += h;
        let mut sm = spec;
        sm.forward_power_price_fp -= h;
        let fd = (evaluate_eydeland_clean_spark_spread_option(&sp).bjerksund_stensland_spark_call
            - evaluate_eydeland_clean_spark_spread_option(&sm).bjerksund_stensland_spark_call)
            / (2.0 * h);
        assert!((r_ad.bjerksund_stensland_spark_call.d[0] - fd).abs() < 1e-5);
    }
}
