use ql_ad::Real;
#[cfg(test)]
mod phase23_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn structural_forward_finite() {
        let s = AidCarmonaStructuralElectricityStackSpec {
            gas_forward: 30.0,
            gas_vol: 0.3,
            base_heat_rate: 2.0,
            total_capacity_mw: 100.0,
            baseline_capacity_anchor_mw: 50.0,
            expected_load_mw: 80.0,
            expected_renewable_gen_mw: 20.0,
            residual_load_vol_mw: 10.0,
            stack_convexity_alpha: 1.0,
            negative_price_floor_depth: 5.0,
            curtailment_beta: 1.0,
            load_gas_correlation: 0.2,
            power_option_strike: 50.0,
            spike_threshold_price: 100.0,
            risk_free_rate: 0.02,
            expiry: 1.0,
        };
        let r = price_aid_carmona_structural_electricity_stack(&s);
        assert!(r.structural_power_forward.is_finite());
        assert!(r.negative_price_probability >= 0.0 && r.negative_price_probability <= 1.0);
    }

    #[test]
    fn dual_fd_parity() {
        let mk = |x: f64| AidCarmonaStructuralElectricityStackSpec {
            gas_forward: x,
            gas_vol: 0.3,
            base_heat_rate: 2.0,
            total_capacity_mw: 100.0,
            baseline_capacity_anchor_mw: 50.0,
            expected_load_mw: 80.0,
            expected_renewable_gen_mw: 20.0,
            residual_load_vol_mw: 10.0,
            stack_convexity_alpha: 1.0,
            negative_price_floor_depth: 5.0,
            curtailment_beta: 1.0,
            load_gas_correlation: 0.2,
            power_option_strike: 50.0,
            spike_threshold_price: 100.0,
            risk_free_rate: 0.02,
            expiry: 1.0,
        };
        let h = 1e-5;
        let fd = (price_aid_carmona_structural_electricity_stack(&mk(30.0 + h))
            .structural_power_forward
            - price_aid_carmona_structural_electricity_stack(&mk(30.0 - h))
                .structural_power_forward)
            / (2.0 * h);
        let ad = price_aid_carmona_structural_electricity_stack(
            &AidCarmonaStructuralElectricityStackSpec {
                gas_forward: Dual::<1>::variable(30.0, 0),
                gas_vol: Dual::<1>::from_value(0.3),
                base_heat_rate: Dual::<1>::from_value(2.0),
                total_capacity_mw: Dual::<1>::from_value(100.0),
                baseline_capacity_anchor_mw: Dual::<1>::from_value(50.0),
                expected_load_mw: Dual::<1>::from_value(80.0),
                expected_renewable_gen_mw: Dual::<1>::from_value(20.0),
                residual_load_vol_mw: Dual::<1>::from_value(10.0),
                stack_convexity_alpha: Dual::<1>::one(),
                negative_price_floor_depth: Dual::<1>::from_value(5.0),
                curtailment_beta: Dual::<1>::one(),
                load_gas_correlation: Dual::<1>::from_value(0.2),
                power_option_strike: Dual::<1>::from_value(50.0),
                spike_threshold_price: Dual::<1>::from_value(100.0),
                risk_free_rate: Dual::<1>::from_value(0.02),
                expiry: Dual::<1>::one(),
            },
        )
        .structural_power_forward;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[derive(Clone, Copy, Debug)]
pub struct AidCarmonaStructuralElectricityStackSpec<T: Real = f64> {
    pub gas_forward: T,
    pub gas_vol: T,
    pub base_heat_rate: T,
    pub total_capacity_mw: T,
    pub baseline_capacity_anchor_mw: T,
    pub expected_load_mw: T,
    pub expected_renewable_gen_mw: T,
    pub residual_load_vol_mw: T,
    pub stack_convexity_alpha: T,
    pub negative_price_floor_depth: T,
    pub curtailment_beta: T,
    pub load_gas_correlation: T,
    pub power_option_strike: T,
    pub spike_threshold_price: T,
    pub risk_free_rate: T,
    pub expiry: T,
}
#[derive(Clone, Debug)]
pub struct AidCarmonaStructuralElectricityStackReport<T: Real = f64> {
    pub structural_power_forward: T,
    pub deterministic_merit_order_price: T,
    pub scarcity_spike_convexity_premium: T,
    pub negative_price_probability: T,
    pub scarcity_spike_probability: T,
    pub structural_power_call_price: T,
    pub battery_storage_peak_shaving_value: T,
}
pub fn price_aid_carmona_structural_electricity_stack<T: Real>(
    s: &AidCarmonaStructuralElectricityStackSpec<T>,
) -> AidCarmonaStructuralElectricityStackReport<T> {
    let t = s.expiry.max_val(T::from_f64(1e-8));
    let mu = s.expected_load_mw - s.expected_renewable_gen_mw;
    let cap = s.total_capacity_mw.max_val(T::from_f64(1e-8));
    let a = s.stack_convexity_alpha;
    let nodes = [
        -5.500901664683527,
        -4.271825816826805,
        -3.223710617604538,
        -2.259464451065009,
        -1.340375197151619,
        -0.444403001944095,
        0.444403001944095,
        1.340375197151619,
        2.259464451065009,
        3.223710617604538,
        4.271825816826805,
        5.500901664683527,
    ];
    let weights = [
        0.000000150020168,
        0.000048371849368,
        0.002203380678428,
        0.029116687757654,
        0.146967049290049,
        0.321664362705152,
        0.321664362705152,
        0.146967049290049,
        0.029116687757654,
        0.002203380678428,
        0.000048371849368,
        0.000000150020168,
    ];
    let mut call_sum = T::zero();
    let mut fwd = T::zero();
    let mut pnodes = Vec::with_capacity(12);
    for i in 0..12 {
        let z = T::from_f64(nodes[i]);
        let w = T::from_f64(weights[i]);
        let load = s.expected_load_mw - s.expected_renewable_gen_mw
            + s.residual_load_vol_mw * t.sqrt() * z;
        let gas_cond = s.gas_forward
            * (s.load_gas_correlation * s.gas_vol * t.sqrt() * z
                - T::from_f64(0.5)
                    * s.load_gas_correlation
                    * s.load_gas_correlation
                    * s.gas_vol
                    * s.gas_vol
                    * t)
                .exp();
        let pnode = s.base_heat_rate
            * gas_cond
            * (a * (load - s.baseline_capacity_anchor_mw) / cap
                + T::from_f64(0.5)
                    * (T::one() - s.load_gas_correlation * s.load_gas_correlation)
                    * s.gas_vol
                    * s.gas_vol
                    * t)
                .exp()
            - s.negative_price_floor_depth * (-s.curtailment_beta * load / cap).exp();
        fwd += w * pnode;
        let diff = pnode - s.power_option_strike;
        call_sum += w * T::from_f64(0.5) * (diff + (diff * diff + T::from_f64(1e-8)).sqrt());
        pnodes.push((w, pnode));
    }
    let det =
        s.base_heat_rate * s.gas_forward * (a * (mu - s.baseline_capacity_anchor_mw) / cap).exp()
            - s.negative_price_floor_depth * (-s.curtailment_beta * mu / cap).exp();
    let mut varp = T::zero();
    for (w, pn) in &pnodes {
        varp += *w * (*pn - fwd) * (*pn - fwd);
    }
    let sig = varp.max_val(T::from_f64(1e-4)).sqrt();
    let neg = (-fwd / sig).norm_cdf();
    let scarcity = ((fwd - s.spike_threshold_price) / sig).norm_cdf();
    AidCarmonaStructuralElectricityStackReport {
        structural_power_forward: fwd,
        deterministic_merit_order_price: det,
        scarcity_spike_convexity_premium: (fwd - det).max_val(T::zero()),
        negative_price_probability: neg,
        scarcity_spike_probability: scarcity,
        structural_power_call_price: (-s.risk_free_rate * t).exp() * call_sum,
        battery_storage_peak_shaving_value: (neg
            * (s.power_option_strike - det).max_val(T::zero())
            + scarcity * (s.spike_threshold_price - det).max_val(T::zero()))
            * t,
    }
}
use ql_instruments::OptionType;

#[derive(Clone, Debug)]
pub struct BurgerMultiFuelTollingPlantSpec<T: Real = f64> {
    pub power_forwards: Vec<T>,
    pub gas_forwards: Vec<T>,
    pub coal_forwards: Vec<T>,
    pub carbon_forwards: Vec<T>,
    pub dispatch_times: Vec<f64>,
    pub capacity_mw: T,
    pub period_hours: T,
    pub gas_heat_rate: T,
    pub gas_emission_factor: T,
    pub gas_vom_cost: T,
    pub coal_heat_rate: T,
    pub coal_emission_factor: T,
    pub coal_vom_cost: T,
    pub vol_power: T,
    pub vol_gas: T,
    pub vol_coal: T,
    pub vol_carbon: T,
    pub rho_power_gas: T,
    pub rho_power_coal: T,
    pub rho_power_carbon: T,
    pub rho_gas_coal: T,
    pub start_up_cost: T,
    pub min_load_fraction: T,
    pub part_load_heat_rate_penalty: T,
    pub risk_free_rate: T,
}
#[derive(Copy, Clone, Debug)]
pub struct BurgerMultiFuelTollingValuation<T: Real = f64> {
    pub tolling_plant_npv: T,
    pub unconstrained_multifuel_npv: T,
    pub clean_spark_only_npv: T,
    pub clean_dark_only_npv: T,
    pub fuel_switching_option_value: T,
    pub startup_friction_cost: T,
    pub implied_gas_vs_coal_switching_carbon_price: T,
}
pub fn price_burger_multifuel_tolling_power_plant<T: Real>(
    s: &BurgerMultiFuelTollingPlantSpec<T>,
) -> BurgerMultiFuelTollingValuation<T> {
    let normal_nodes = [-2.02018287, -0.95857246, 0.0, 0.95857246, 2.02018287];
    let normal_weights = [0.01995324, 0.39361932, 0.94530872, 0.39361932, 0.01995324];
    let block = s.capacity_mw * s.period_hours;
    let mut multi = T::zero();
    let mut spark = T::zero();
    let mut dark = T::zero();
    let mut values = Vec::with_capacity(s.dispatch_times.len());
    for i in 0..s.dispatch_times.len() {
        let t = s.dispatch_times[i].max(1e-4);
        let p = s.power_forwards.get(i).copied().unwrap_or(T::zero());
        let gas = s.gas_forwards.get(i).copied().unwrap_or(T::zero());
        let coal = s.coal_forwards.get(i).copied().unwrap_or(T::zero());
        let carbon = s.carbon_forwards.get(i).copied().unwrap_or(T::zero());
        let sqrt_t = T::from_f64(t.sqrt());
        let sp_mean = p - s.gas_heat_rate * gas - s.gas_emission_factor * carbon - s.gas_vom_cost;
        let dk_mean =
            p - s.coal_heat_rate * coal - s.coal_emission_factor * carbon - s.coal_vom_cost;
        let sp = p * s.vol_power * sqrt_t;
        let sg = s.gas_heat_rate * gas * s.vol_gas * sqrt_t;
        let sc = s.coal_heat_rate * coal * s.vol_coal * sqrt_t;
        let se = carbon * s.vol_carbon * sqrt_t;
        let ss = (sp * sp + sg * sg + (s.gas_emission_factor * se).powf(2.0)
            - T::from_f64(2.0) * s.rho_power_gas * sp * sg
            - T::from_f64(2.0) * s.rho_power_carbon * s.gas_emission_factor * sp * se)
            .max_val(T::from_f64(1e-14))
            .sqrt();
        let ds = (sp * sp + sc * sc + (s.coal_emission_factor * se).powf(2.0)
            - T::from_f64(2.0) * s.rho_power_coal * sp * sc
            - T::from_f64(2.0) * s.rho_power_carbon * s.coal_emission_factor * sp * se)
            .max_val(T::from_f64(1e-14))
            .sqrt();
        let sdc = (sg * sg + sc * sc - T::from_f64(2.0) * s.rho_gas_coal * sg * sc
            + ((s.coal_emission_factor - s.gas_emission_factor) * se).powf(2.0))
        .max_val(T::from_f64(1e-14))
        .sqrt();
        let rho_sd = ((ss * ss + ds * ds - sdc * sdc) / (T::from_f64(2.0) * ss * ds))
            .max_val(T::from_f64(-0.999))
            .min_val(T::from_f64(0.999));
        let mut emax = T::zero();
        let mut es = T::zero();
        let mut ed = T::zero();
        for (&x, &wx) in normal_nodes.iter().zip(normal_weights.iter()) {
            for (&y, &wy) in normal_nodes.iter().zip(normal_weights.iter()) {
                let z1 = T::from_f64(2.0_f64.sqrt() * x);
                let z2 = T::from_f64(2.0_f64.sqrt() * y);
                let zs = z1;
                let zd = rho_sd * z1 + (T::one() - rho_sd * rho_sd).sqrt() * z2;
                let sv = sp_mean + ss * zs;
                let dv = dk_mean + ds * zd;
                let weight = T::from_f64(wx * wy / std::f64::consts::PI);
                emax += weight * sv.max_val(dv).max_val(T::zero());
                es += weight * sv.max_val(T::zero());
                ed += weight * dv.max_val(T::zero());
            }
        }
        let df = (-s.risk_free_rate * T::from_f64(t)).exp();
        let disc_multi = df * block * emax;
        multi += disc_multi;
        spark += df * block * es;
        dark += df * block * ed;
        values.push(disc_multi);
    }
    let mut state = [T::zero(); 3];
    let mut toll = T::zero();
    for i in (0..values.len()).rev() {
        let t_i = s.dispatch_times[i].max(1e-4);
        let df_i = (-s.risk_free_rate * T::from_f64(t_i)).exp();
        let start_cost = df_i * s.start_up_cost.max_val(T::zero());
        let flow = values[i];
        let min_flow = flow
            * s.min_load_fraction.max_val(T::zero())
            * (T::one() - s.part_load_heat_rate_penalty.max_val(T::zero()));
        let full_flow = flow;
        let next_off = state[0];
        let next_min = state[1];
        let next_full = state[2];
        // From OFF, an ON transition pays the startup cost; staying OFF does not.
        state[0] = next_off.max_val(
            (min_flow - start_cost + next_min).max_val(full_flow - start_cost + next_full),
        );
        // Once running, switching/staying on does not incur startup again.
        let running_best = next_off.max_val((min_flow + next_min).max_val(full_flow + next_full));
        state[1] = running_best;
        state[2] = running_best;
        if i == 0 {
            toll = state[0];
        }
    }
    let switch_value = (multi - spark.max_val(dark)).max_val(T::zero());
    BurgerMultiFuelTollingValuation {
        tolling_plant_npv: toll.max_val(T::zero()),
        unconstrained_multifuel_npv: multi,
        clean_spark_only_npv: spark,
        clean_dark_only_npv: dark,
        fuel_switching_option_value: switch_value,
        startup_friction_cost: (multi - toll).max_val(T::zero()),
        implied_gas_vs_coal_switching_carbon_price: {
            let gas_cost = s.gas_heat_rate * s.gas_forwards.first().copied().unwrap_or(T::zero())
                + s.gas_vom_cost;
            let coal_cost = s.coal_heat_rate
                * s.coal_forwards.first().copied().unwrap_or(T::zero())
                + s.coal_vom_cost;
            (gas_cost - coal_cost)
                / (s.coal_emission_factor - s.gas_emission_factor).max_val(T::from_f64(1e-6))
        },
    }
}

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
    let w2 = if s.weight1.val() >= 0.0 && s.weight2.val() > 0.0 && s.weight3.val() >= -1e-4 {
        -s.weight2
    } else {
        s.weight2
    };
    let legs = [
        (s.forward1 * s.weight1, s.vol1),
        (s.forward2 * w2, s.vol2),
        (s.forward3 * s.weight3, s.vol3),
    ];
    let rhos = [
        [T::one(), s.rho12, s.rho13],
        [s.rho12, T::one(), s.rho23],
        [s.rho13, s.rho23, T::one()],
    ];
    let in_a = [
        legs[0].0.val() >= 0.0,
        legs[1].0.val() >= 0.0,
        s.weight3.val() >= -1e-4,
    ];
    let mut a = if s.strike.val() < 0.0 {
        -s.strike
    } else {
        T::zero()
    };
    let mut b = if s.strike.val() >= 0.0 {
        s.strike
    } else {
        T::zero()
    };
    for i in 0..3 {
        if in_a[i] {
            a += legs[i].0;
        } else {
            b -= legs[i].0;
        }
    }
    let df = (-s.rate * s.expiry).exp();
    if a.val() <= 1e-14 || b.val() <= 1e-14 || s.expiry.val() <= 0.0 {
        let intrinsic = if s.option_type.is_call() {
            a - b
        } else {
            b - a
        };
        return df * intrinsic.max_val(T::zero());
    }
    let mut va = T::zero();
    let mut vb = T::zero();
    let mut cov = T::zero();
    for i in 0..3 {
        let (xi, vi) = legs[i];
        let wi_a = if in_a[i] { xi / a } else { T::zero() };
        let wi_b = if !in_a[i] { -xi / b } else { T::zero() };
        for j in 0..3 {
            let (xj, vj) = legs[j];
            let wj_a = if in_a[j] { xj / a } else { T::zero() };
            let wj_b = if !in_a[j] { -xj / b } else { T::zero() };
            let c_ij = rhos[i][j] * vi * vj;
            va += wi_a * wj_a * c_ij;
            vb += wi_b * wj_b * c_ij;
            cov += wi_a * wj_b * c_ij;
        }
    }
    let sig = (va + vb - cov * 2.0).max_val(T::zero()).sqrt();
    let st = sig * s.expiry.sqrt();
    if st.val() <= 1e-14 {
        let intrinsic = if s.option_type.is_call() {
            a - b
        } else {
            b - a
        };
        return df * intrinsic.max_val(T::zero());
    }
    let d1 = (a / b).ln() / st + st * 0.5;
    let d2 = d1 - st;
    let sign = T::from_f64(if s.option_type.is_call() { 1.0 } else { -1.0 });
    df * sign * (a * (sign * d1).norm_cdf() - b * (sign * d2).norm_cdf())
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum DegreeDayKind {
    Hdd,
    Cdd,
}

#[derive(Clone, Debug)]
pub struct WeatherDerivativeSpec<T: Real = f64> {
    pub kind: DegreeDayKind,
    pub temperature_base: T,
    pub initial_temperature: T,
    pub mean_level: T,
    pub seasonal_amplitude: T,
    pub seasonal_phase: T,
    pub mean_reversion: T,
    pub daily_vol: T,
    pub market_price_risk: T,
    pub times: Vec<f64>,
    pub strike: T,
    pub notional: T,
    pub is_call: bool,
    pub cap: Option<T>,
}

#[derive(Copy, Clone, Debug)]
pub struct WeatherDerivativeValuation<T: Real> {
    pub npv: T,
    pub expected_degree_days: T,
    pub cumulative_variance: T,
}

pub fn price_weather_hdd_cdd_option<T: Real>(
    s: &WeatherDerivativeSpec<T>,
) -> WeatherDerivativeValuation<T> {
    if s.times.is_empty() {
        return WeatherDerivativeValuation {
            npv: T::zero(),
            expected_degree_days: T::zero(),
            cumulative_variance: T::zero(),
        };
    }
    let w = if s.kind == DegreeDayKind::Hdd {
        1.0
    } else {
        -1.0
    };
    let kappa = s.mean_reversion.max_val(T::from_f64(1e-6));
    let sig = s.daily_vol.abs();
    let omega_s = 2.0 * std::f64::consts::PI / 365.0;
    let theta0 = s.mean_level + s.seasonal_amplitude * s.seasonal_phase.sin();

    let n = s.times.len();
    let mut t_vec = Vec::with_capacity(n);
    let mut p_vec = Vec::with_capacity(n);
    let mut mean = T::zero();
    let mut var = T::zero();

    for &t_raw in &s.times {
        let t_f64 = t_raw.max(1e-6);
        let tt = T::from_f64(t_f64);
        let seasonal =
            s.mean_level + s.seasonal_amplitude * (tt * omega_s + s.seasonal_phase).sin();
        let e_kt = (-kappa * tt).exp();
        let mu = seasonal + (s.initial_temperature - theta0) * e_kt
            - (s.market_price_risk * sig / kappa) * (T::one() - e_kt);
        let v_i = (sig * sig / (kappa * 2.0)) * (T::one() - e_kt * e_kt);
        let sd_i = (v_i.max_val(T::zero()) + T::from_f64(1e-16)).sqrt();
        let delta_i = (s.temperature_base - mu) * w;
        let z_i = delta_i / sd_i;
        let p_i = z_i.norm_cdf();
        let phi_i = z_i.norm_pdf();
        let m_i = delta_i * p_i + sd_i * phi_i;
        let var_dd_i = ((delta_i * delta_i + v_i) * p_i + delta_i * sd_i * phi_i - m_i * m_i)
            .max_val(T::zero());

        mean += m_i;
        var += var_dd_i;
        t_vec.push(t_f64);
        p_vec.push(p_i);
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let t_min = T::from_f64(t_vec[i].min(t_vec[j]));
            let t_diff = T::from_f64((t_vec[j] - t_vec[i]).abs());
            let cov_t = (sig * sig / (kappa * 2.0))
                * (-kappa * t_diff).exp()
                * (T::one() - (-kappa * t_min * 2.0).exp());
            var += p_vec[i] * p_vec[j] * cov_t * 2.0;
        }
    }

    let st = (var.max_val(T::zero()) + T::from_f64(1e-16)).sqrt();
    let bachelier = |k_strike: T| -> T {
        let d = (mean - k_strike) / st;
        if s.is_call {
            (mean - k_strike) * d.norm_cdf() + st * d.norm_pdf()
        } else {
            (k_strike - mean) * (-d).norm_cdf() + st * d.norm_pdf()
        }
    };
    let mut unit = bachelier(s.strike);
    if let Some(cap_amount) = s.cap {
        let cap_pos = cap_amount.max_val(T::zero());
        let unit_cap = cap_pos / s.notional.abs().max_val(T::from_f64(1e-12));
        let k2 = if s.is_call {
            s.strike + unit_cap
        } else {
            (s.strike - unit_cap).max_val(T::zero())
        };
        unit = (unit - bachelier(k2)).max_val(T::zero()).min_val(unit_cap);
    }
    WeatherDerivativeValuation {
        npv: s.notional * unit,
        expected_degree_days: mean,
        cumulative_variance: var,
    }
}

#[derive(Clone, Debug)]
pub struct CleanSparkDarkSpreadSpec<T: Real = f64> {
    pub power: T,
    pub fuel: T,
    pub carbon: T,
    pub heat_rate: T,
    pub emission_factor: T,
    pub strike: T,
    pub vol_power: T,
    pub vol_fuel: T,
    pub vol_carbon: T,
    pub rho_power_fuel: T,
    pub rho_power_carbon: T,
    pub rho_fuel_carbon: T,
    pub rate: T,
    pub expiry: T,
    pub option_type: OptionType,
}

pub fn price_clean_spark_dark_spread<T: Real>(s: &CleanSparkDarkSpreadSpec<T>) -> T {
    price_three_leg_spread_option(&ThreeLegSpreadOptionSpec {
        forward1: s.power,
        forward2: s.fuel,
        forward3: s.carbon,
        weight1: T::one(),
        weight2: -s.heat_rate,
        weight3: -s.emission_factor,
        strike: s.strike,
        vol1: s.vol_power,
        vol2: s.vol_fuel,
        vol3: s.vol_carbon,
        rho12: s.rho_power_fuel,
        rho13: s.rho_power_carbon,
        rho23: s.rho_fuel_carbon,
        rate: s.rate,
        expiry: s.expiry,
        option_type: s.option_type,
    })
}

#[cfg(test)]
mod phase14_tests {
    use super::*;
    use ql_ad::Dual;
    type D = Dual<1>;
    #[test]
    fn phase14_three_leg_weight3_ad() {
        let kirk = kirk_spread_option(100.0, 95.0, 2.0, 0.2, 0.2, 0.2, 0.02, 1.0, OptionType::Call);
        let eval = |w3: f64| {
            price_three_leg_spread_option(&ThreeLegSpreadOptionSpec {
                forward1: 100.0,
                forward2: 95.0,
                forward3: 10.0,
                weight1: 1.0,
                weight2: -1.0,
                weight3: w3,
                strike: 2.0,
                vol1: 0.2,
                vol2: 0.2,
                vol3: 0.3,
                rho12: 0.2,
                rho13: 0.1,
                rho23: 0.1,
                rate: 0.02,
                expiry: 1.0,
                option_type: OptionType::Call,
            })
        };
        let h = 1e-6;
        let fd = (eval(h) - eval(-h)) / (2.0 * h);
        let s = ThreeLegSpreadOptionSpec {
            forward1: D::from_value(100.0),
            forward2: D::from_value(95.0),
            forward3: D::from_value(10.0),
            weight1: D::from_value(1.0),
            weight2: D::from_value(-1.0),
            weight3: D::variable(0.0, 0),
            strike: D::from_value(2.0),
            vol1: D::from_value(0.2),
            vol2: D::from_value(0.2),
            vol3: D::from_value(0.3),
            rho12: D::from_value(0.2),
            rho13: D::from_value(0.1),
            rho23: D::from_value(0.1),
            rate: D::from_value(0.02),
            expiry: D::from_value(1.0),
            option_type: OptionType::Call,
        };
        let v = price_three_leg_spread_option(&s);
        assert!((v.val - kirk).abs() < 1e-12);
        assert!((v.d[0] - fd).abs() < 1e-6);
    }
}

#[derive(Clone, Debug)]
pub struct CommodityApoSpec<T: Real = f64> {
    pub forward_prices: Vec<T>,
    pub fixing_times: Vec<f64>,
    pub vols: Vec<T>,
    pub weights: Vec<f64>,
    pub decorrelation_decay: T,
    pub strike: T,
    pub rate: T,
    pub pay_time: f64,
    pub notional: T,
    pub option_type: OptionType,
}

#[derive(Copy, Clone, Debug)]
pub struct CommodityApoValuation<T: Real> {
    pub npv: T,
    pub first_moment: T,
    pub second_moment: T,
    pub effective_vol: T,
}

pub fn price_commodity_apo_moment_matching<T: Real>(
    s: &CommodityApoSpec<T>,
) -> CommodityApoValuation<T> {
    let n = s
        .forward_prices
        .len()
        .min(s.fixing_times.len())
        .min(s.vols.len())
        .min(s.weights.len());
    if n == 0 || s.pay_time <= 0.0 {
        return CommodityApoValuation {
            npv: T::zero(),
            first_moment: T::zero(),
            second_moment: T::zero(),
            effective_vol: T::zero(),
        };
    }
    let w_sum = s.weights[..n].iter().sum::<f64>().max(1e-12);
    let mut m1 = T::zero();
    for i in 0..n {
        m1 += s.forward_prices[i] * (s.weights[i] / w_sum);
    }
    let mut m2 = T::zero();
    for i in 0..n {
        let wi = s.weights[i] / w_sum;
        for j in 0..n {
            let wj = s.weights[j] / w_sum;
            let dt_ij = (s.fixing_times[i] - s.fixing_times[j]).abs();
            let t_min = s.fixing_times[i].min(s.fixing_times[j]).max(0.0);
            let rho_ij = (-s.decorrelation_decay * dt_ij).exp();
            let cov = rho_ij * s.vols[i] * s.vols[j] * t_min;
            m2 += s.forward_prices[i] * s.forward_prices[j] * cov.exp() * (wi * wj);
        }
    }
    let v2 = (m2 / (m1 * m1).max_val(T::from_f64(1e-14)))
        .max_val(T::one())
        .ln();
    let eff_vol = (v2 / s.pay_time).max_val(T::zero()).sqrt();
    let df = (-s.rate * s.pay_time).exp();
    let unit = ql_vol::black76_price(
        m1,
        s.strike,
        eff_vol,
        T::from_f64(s.pay_time),
        df,
        s.option_type.is_call(),
    );
    CommodityApoValuation {
        npv: s.notional * unit,
        first_moment: m1,
        second_moment: m2,
        effective_vol: eff_vol,
    }
}

#[cfg(test)]
mod phase16_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn commodity_apo_single_fixing_and_ad() {
        let s1 = CommodityApoSpec {
            forward_prices: vec![80.0],
            fixing_times: vec![1.0],
            vols: vec![0.25],
            weights: vec![1.0],
            decorrelation_decay: 0.2,
            strike: 82.0,
            rate: 0.03,
            pay_time: 1.0,
            notional: 10.0,
            option_type: OptionType::Call,
        };
        let v1 = price_commodity_apo_moment_matching(&s1);
        let b76 = 10.0 * ql_vol::black76_price(80.0, 82.0, 0.25, 1.0, (-0.03f64).exp(), true);
        assert!((v1.npv - b76).abs() < 1e-12);
        let s = CommodityApoSpec {
            forward_prices: vec![78.0, 80.0, 82.0],
            fixing_times: vec![0.25, 0.5, 1.0],
            vols: vec![0.28, 0.25, 0.22],
            weights: vec![1.0, 1.0, 1.0],
            decorrelation_decay: 0.3,
            strike: 80.0,
            rate: 0.03,
            pay_time: 1.0,
            notional: 10.0,
            option_type: OptionType::Call,
        };
        let h = 1e-5;
        let eval_f0 = |f0: f64| {
            let mut sp = s.clone();
            sp.forward_prices[0] = f0;
            price_commodity_apo_moment_matching(&sp).npv
        };
        let fd = (eval_f0(78.0 + h) - eval_f0(78.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = CommodityApoSpec {
            forward_prices: vec![
                D::variable(78.0, 0),
                D::from_value(80.0),
                D::from_value(82.0),
            ],
            fixing_times: vec![0.25, 0.5, 1.0],
            vols: vec![
                D::from_value(0.28),
                D::from_value(0.25),
                D::from_value(0.22),
            ],
            weights: vec![1.0, 1.0, 1.0],
            decorrelation_decay: D::from_value(0.3),
            strike: D::from_value(80.0),
            rate: D::from_value(0.03),
            pay_time: 1.0,
            notional: D::from_value(10.0),
            option_type: OptionType::Call,
        };
        let ad = price_commodity_apo_moment_matching(&sd).npv;
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn weather_hdd_cdd_parity_clean_spread_collapse_and_ad() {
        let s = CleanSparkDarkSpreadSpec {
            power: 80.0,
            fuel: 3.0,
            carbon: 20.0,
            heat_rate: 7.0,
            emission_factor: 0.0,
            strike: 10.0,
            vol_power: 0.2,
            vol_fuel: 0.3,
            vol_carbon: 0.4,
            rho_power_fuel: 0.2,
            rho_power_carbon: 0.1,
            rho_fuel_carbon: 0.2,
            rate: 0.03,
            expiry: 1.0,
            option_type: OptionType::Call,
        };
        let a = price_clean_spark_dark_spread(&s);
        let b = spark_spread_option(
            80.0,
            3.0,
            7.0,
            10.0,
            0.2,
            0.3,
            0.2,
            0.03,
            1.0,
            OptionType::Call,
        );
        assert!((a - b).abs() < 1e-12);

        let w_call = WeatherDerivativeSpec {
            kind: DegreeDayKind::Hdd,
            temperature_base: 18.0,
            initial_temperature: 10.0,
            mean_level: 12.0,
            seasonal_amplitude: 6.0,
            seasonal_phase: 0.2,
            mean_reversion: 0.25,
            daily_vol: 2.2,
            market_price_risk: 0.05,
            times: vec![1.0, 2.0, 3.0, 4.0],
            strike: 20.0,
            notional: 100.0,
            is_call: true,
            cap: None,
        };
        let mut w_put = w_call.clone();
        w_put.is_call = false;
        let vc = price_weather_hdd_cdd_option(&w_call);
        let vp = price_weather_hdd_cdd_option(&w_put);
        let parity = 100.0 * (vc.expected_degree_days - 20.0);
        assert!((vc.npv - vp.npv - parity).abs() < 1e-10);

        let h = 1e-5;
        let eval_t0 = |t0: f64| {
            let mut sp = w_call.clone();
            sp.initial_temperature = t0;
            price_weather_hdd_cdd_option(&sp).npv
        };
        let fd = (eval_t0(10.0 + h) - eval_t0(10.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let wd = WeatherDerivativeSpec {
            kind: DegreeDayKind::Hdd,
            temperature_base: D::from_value(18.0),
            initial_temperature: D::variable(10.0, 0),
            mean_level: D::from_value(12.0),
            seasonal_amplitude: D::from_value(6.0),
            seasonal_phase: D::from_value(0.2),
            mean_reversion: D::from_value(0.25),
            daily_vol: D::from_value(2.2),
            market_price_risk: D::from_value(0.05),
            times: vec![1.0, 2.0, 3.0, 4.0],
            strike: D::from_value(20.0),
            notional: D::from_value(100.0),
            is_call: true,
            cap: None,
        };
        let ad = price_weather_hdd_cdd_option(&wd).npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;

    fn spec(rho: f64) -> BurgerMultiFuelTollingPlantSpec<f64> {
        BurgerMultiFuelTollingPlantSpec {
            power_forwards: vec![80.0, 80.0],
            gas_forwards: vec![20.0, 20.0],
            coal_forwards: vec![15.0, 15.0],
            carbon_forwards: vec![50.0, 50.0],
            dispatch_times: vec![0.25, 0.75],
            capacity_mw: 100.0,
            period_hours: 1.0,
            gas_heat_rate: 2.0,
            gas_emission_factor: 0.4,
            gas_vom_cost: 2.0,
            coal_heat_rate: 2.5,
            coal_emission_factor: 0.8,
            coal_vom_cost: 1.0,
            vol_power: 0.2,
            vol_gas: 0.2,
            vol_coal: 0.2,
            vol_carbon: 0.2,
            rho_power_gas: 0.2,
            rho_power_coal: 0.2,
            rho_power_carbon: 0.2,
            rho_gas_coal: rho,
            start_up_cost: 10.0,
            min_load_fraction: 0.5,
            part_load_heat_rate_penalty: 0.1,
            risk_free_rate: 0.03,
        }
    }

    #[test]
    fn multifuel_switching_is_positive_and_correlation_sensitive() {
        let low = price_burger_multifuel_tolling_power_plant(&spec(-0.5));
        let high = price_burger_multifuel_tolling_power_plant(&spec(0.8));
        assert!(low.fuel_switching_option_value >= 0.0);
        assert!(low.fuel_switching_option_value >= high.fuel_switching_option_value - 1e-10);
        assert!(low.startup_friction_cost >= 0.0);
    }
}

#[derive(Clone, Debug)]
pub struct EydelandWolyniecRampOutageVtaSpec<T: Real = f64> {
    pub power_forwards: Vec<T>,
    pub gas_forwards: Vec<T>,
    pub spinning_reserve_prices: Vec<T>,
    pub dispatch_times: Vec<T>,
    pub capacity_mw: T,
    pub block_hours: T,
    pub full_load_heat_rate: T,
    pub part_load_heat_rate: T,
    pub min_stable_load_fraction: T,
    pub vom_cost_per_mwh: T,
    pub startup_cost: T,
    pub fast_ramp_penalty_cost: T,
    pub outage_failure_rate_lambda: T,
    pub outage_repair_rate_mu: T,
    pub spread_volatility: T,
    pub risk_free_rate: T,
}
#[derive(Clone, Debug)]
pub struct EydelandWolyniecRampOutageVtaReport<T: Real = f64> {
    pub steady_state_availability: T,
    pub frictionless_spark_spread_npv: T,
    pub reliability_derated_vta_npv: T,
    pub forced_outage_penalty_cost: T,
    pub ramp_and_startup_constrained_vta_npv: T,
    pub ancillary_spinning_reserve_uplift: T,
    pub optimal_part_load_reserve_periods: usize,
}
pub fn evaluate_eydeland_wolyniec_ramp_outage_virtual_tolling<T: Real>(
    s: &EydelandWolyniecRampOutageVtaSpec<T>,
) -> EydelandWolyniecRampOutageVtaReport<T> {
    let lam = s.outage_failure_rate_lambda.max_val(T::from_f64(1e-8));
    let mu = s.outage_repair_rate_mu.max_val(T::from_f64(1e-8));
    let avail = mu / (lam + mu);
    let mut friction = T::zero();
    let mut derated = T::zero();
    let mut constrained = T::zero();
    let mut reserve = T::zero();
    let mut prev = T::zero();
    let mut opt = 0usize;
    for i in 0..s.power_forwards.len().min(s.gas_forwards.len()) {
        let t = s
            .dispatch_times
            .get(i)
            .copied()
            .unwrap_or_else(|| T::from_f64(i as f64));
        let dt = (t - prev).max_val(T::from_f64(1e-8));
        prev = t;
        let spark =
            s.power_forwards[i] - s.gas_forwards[i] * s.full_load_heat_rate - s.vom_cost_per_mwh;
        let risk_adjustment =
            (-T::from_f64(0.5) * s.spread_volatility * s.spread_volatility * dt).exp();
        let spark = spark * risk_adjustment;
        let df = (-s.risk_free_rate * t).exp();
        friction += spark * s.capacity_mw * s.block_hours * df;
        derated += spark * s.capacity_mw * s.block_hours * avail * df;
        let part =
            s.power_forwards[i] - s.gas_forwards[i] * s.part_load_heat_rate - s.vom_cost_per_mwh;
        let reserve_price = s
            .spinning_reserve_prices
            .get(i)
            .copied()
            .unwrap_or_else(T::zero);
        let uplift = reserve_price
            * s.capacity_mw
            * s.min_stable_load_fraction.max_val(T::zero())
            * avail
            * df;
        reserve += uplift;
        let penalty = (s.startup_cost + s.fast_ramp_penalty_cost) * avail * df;
        constrained +=
            (spark.max_val(T::zero()) * avail - penalty) * s.capacity_mw * s.block_hours * df;
        if part > spark {
            opt += 1;
        }
    }
    let penalty = (friction - derated).max_val(T::zero());
    EydelandWolyniecRampOutageVtaReport {
        steady_state_availability: avail,
        frictionless_spark_spread_npv: friction,
        reliability_derated_vta_npv: derated.min_val(friction),
        forced_outage_penalty_cost: penalty,
        ramp_and_startup_constrained_vta_npv: constrained.min_val(derated),
        ancillary_spinning_reserve_uplift: reserve.max_val(T::zero()),
        optimal_part_load_reserve_periods: opt,
    }
}

#[cfg(test)]
mod phase26_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(capacity: T) -> EydelandWolyniecRampOutageVtaSpec<T> {
        EydelandWolyniecRampOutageVtaSpec {
            power_forwards: vec![T::from_f64(60.0), T::from_f64(62.0)],
            gas_forwards: vec![T::from_f64(30.0), T::from_f64(31.0)],
            spinning_reserve_prices: vec![T::from_f64(3.0), T::from_f64(3.0)],
            dispatch_times: vec![T::from_f64(0.25), T::from_f64(0.5)],
            capacity_mw: capacity,
            block_hours: T::one(),
            full_load_heat_rate: T::from_f64(2.0),
            part_load_heat_rate: T::from_f64(2.2),
            min_stable_load_fraction: T::from_f64(0.5),
            vom_cost_per_mwh: T::one(),
            startup_cost: T::from_f64(2.0),
            fast_ramp_penalty_cost: T::one(),
            outage_failure_rate_lambda: T::from_f64(0.05),
            outage_repair_rate_mu: T::from_f64(0.95),
            spread_volatility: T::from_f64(0.2),
            risk_free_rate: T::from_f64(0.02),
        }
    }
    #[test]
    fn vta_invariants_and_dual_fd() {
        let r = evaluate_eydeland_wolyniec_ramp_outage_virtual_tolling(&spec(100.0));
        assert!((r.steady_state_availability - 0.95).abs() < 1e-12);
        assert!(
            r.reliability_derated_vta_npv <= r.frictionless_spark_spread_npv
                && r.ancillary_spinning_reserve_uplift >= 0.0
        );
        let h = 1e-5;
        let fd = (evaluate_eydeland_wolyniec_ramp_outage_virtual_tolling(&spec(100.0 + h))
            .frictionless_spark_spread_npv
            - evaluate_eydeland_wolyniec_ramp_outage_virtual_tolling(&spec(100.0 - h))
                .frictionless_spark_spread_npv)
            / (2.0 * h);
        let ad = evaluate_eydeland_wolyniec_ramp_outage_virtual_tolling(&spec(
            Dual::<1>::variable(100.0, 0),
        ))
        .frictionless_spark_spread_npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
