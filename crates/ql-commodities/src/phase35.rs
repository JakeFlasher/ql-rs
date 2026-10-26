use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Clone, Copy, Debug)]
pub struct EydelandGibsonSchwartzKirkSpreadSpec<T: Real = f64> {
    pub spot_price: T,
    pub initial_convenience_yield: T,
    pub mean_reversion_kappa: T,
    pub long_run_convenience_yield: T,
    pub spot_volatility: T,
    pub convenience_yield_vol: T,
    pub spot_cy_correlation: T,
    pub risk_free_rate: T,
    pub maturity_years: T,
    pub fuel_futures_price: T,
    pub heat_rate_multiplier: T,
    pub fixed_strike_cost: T,
    pub fuel_volatility: T,
    pub power_fuel_correlation: T,
}

#[derive(Clone, Copy, Debug)]
pub struct EydelandGibsonSchwartzKirkSpreadReport<T: Real = f64> {
    pub convenience_yield_loading_b: T,
    pub gibson_schwartz_futures_price: T,
    pub integrated_futures_variance: T,
    pub terminal_samuelson_futures_vol: T,
    pub kirk_composite_spread_vol: T,
    pub kirk_spark_spread_call_pv: T,
    pub intrinsic_spark_spread_pv: T,
}

pub fn evaluate_eydeland_gibson_schwartz_and_kirk_spread<T: Real>(
    s: &EydelandGibsonSchwartzKirkSpreadSpec<T>,
) -> EydelandGibsonSchwartzKirkSpreadReport<T> {
    let t = s.maturity_years.max_val(T::from_f64(1e-6));
    let kappa = s.mean_reversion_kappa.max_val(T::from_f64(1e-6));
    let e_kt = (-kappa * t).exp();
    let b = (T::one() - e_kt) / kappa;
    let ss2 = s.spot_volatility * s.spot_volatility;
    let sd2 = s.convenience_yield_vol * s.convenience_yield_vol;
    let v_gs = ss2 * t
        - (T::from_f64(2.0) * s.spot_cy_correlation * s.spot_volatility * s.convenience_yield_vol
            / kappa)
            * (t - b)
        + (sd2 / (T::from_f64(2.0) * kappa * kappa * kappa))
            * (T::from_f64(2.0) * kappa * t - T::from_f64(3.0) + T::from_f64(4.0) * e_kt
                - e_kt * e_kt);
    let f_gs = s.spot_price
        * (-s.initial_convenience_yield * b
            + (s.risk_free_rate - s.long_run_convenience_yield) * t
            + s.long_run_convenience_yield * b
            + T::from_f64(0.5) * (v_gs - ss2 * t))
            .exp();
    let sigma_f = (ss2
        - T::from_f64(2.0)
            * s.spot_cy_correlation
            * s.spot_volatility
            * s.convenience_yield_vol
            * b
        + sd2 * b * b)
        .max_val(T::from_f64(1e-8))
        .sqrt();
    let f_leg2 = s.heat_rate_multiplier * s.fuel_futures_price + s.fixed_strike_cost;
    let w_fuel =
        (s.heat_rate_multiplier * s.fuel_futures_price) / f_leg2.max_val(T::from_f64(1e-8));
    let sigma_leg2 = w_fuel * s.fuel_volatility;
    let sigma_kirk = (sigma_f * sigma_f
        - T::from_f64(2.0) * s.power_fuel_correlation * sigma_f * sigma_leg2
        + sigma_leg2 * sigma_leg2)
        .max_val(T::from_f64(1e-8))
        .sqrt();
    let st = sigma_kirk * t.sqrt();
    let d1 = ((f_gs / f_leg2).ln() + T::from_f64(0.5) * sigma_kirk * sigma_kirk * t) / st;
    let d2 = d1 - st;
    let df = (-s.risk_free_rate * t).exp();
    let kirk_call = df * (f_gs * norm_cdf(d1) - f_leg2 * norm_cdf(d2));
    let diff = f_gs - f_leg2;
    let intr = df * T::from_f64(0.5) * (diff + (diff * diff + T::from_f64(1e-16)).sqrt());
    EydelandGibsonSchwartzKirkSpreadReport {
        convenience_yield_loading_b: b,
        gibson_schwartz_futures_price: f_gs,
        integrated_futures_variance: v_gs,
        terminal_samuelson_futures_vol: sigma_f,
        kirk_composite_spread_vol: sigma_kirk,
        kirk_spark_spread_call_pv: kirk_call,
        intrinsic_spark_spread_pv: intr,
    }
}

#[cfg(test)]
mod phase35_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> EydelandGibsonSchwartzKirkSpreadSpec<T> {
        EydelandGibsonSchwartzKirkSpreadSpec {
            spot_price: T::from_f64(55.0) + x,
            initial_convenience_yield: T::from_f64(0.06),
            mean_reversion_kappa: T::from_f64(1.2),
            long_run_convenience_yield: T::from_f64(0.05),
            spot_volatility: T::from_f64(0.32),
            convenience_yield_vol: T::from_f64(0.20),
            spot_cy_correlation: T::from_f64(0.4),
            risk_free_rate: T::from_f64(0.03),
            maturity_years: T::from_f64(1.0),
            fuel_futures_price: T::from_f64(6.0),
            heat_rate_multiplier: T::from_f64(7.5),
            fixed_strike_cost: T::from_f64(5.0),
            fuel_volatility: T::from_f64(0.28),
            power_fuel_correlation: T::from_f64(0.6),
        }
    }
    #[test]
    fn phase35_commodities_invariants_and_dual_fd() {
        let r = evaluate_eydeland_gibson_schwartz_and_kirk_spread(&spec(0.0));
        assert!(r.kirk_spark_spread_call_pv >= r.intrinsic_spark_spread_pv - 1e-6);
        assert!(r.integrated_futures_variance > 0.0);
        let h = 1e-6;
        let fd = (evaluate_eydeland_gibson_schwartz_and_kirk_spread(&spec(h))
            .kirk_spark_spread_call_pv
            - evaluate_eydeland_gibson_schwartz_and_kirk_spread(&spec(-h))
                .kirk_spark_spread_call_pv)
            / (2.0 * h);
        let ad =
            evaluate_eydeland_gibson_schwartz_and_kirk_spread(&spec(Dual::<1>::variable(0.0, 0)))
                .kirk_spark_spread_call_pv;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
