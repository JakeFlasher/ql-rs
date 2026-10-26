use ql_ad::Real;
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
    let n = s.power_forwards.len().min(s.gas_forwards.len());
    let mut friction = T::zero();
    let mut derated = T::zero();
    let mut reserve = T::zero();
    let mut opt = 0usize;
    let mut v_dp = [T::zero(); 3];
    for i in (0..n).rev() {
        let t = s
            .dispatch_times
            .get(i)
            .copied()
            .unwrap_or_else(|| T::from_f64((i + 1) as f64 * 0.25))
            .max_val(T::from_f64(1e-8));
        let p_avail_t = avail + (lam / (lam + mu)) * (-(lam + mu) * t).exp();
        let df = (-s.risk_free_rate * t).exp();
        let mwh = s.capacity_mw * s.block_hours;
        let s_full =
            s.power_forwards[i] - s.gas_forwards[i] * s.full_load_heat_rate - s.vom_cost_per_mwh;
        let sig_t = (s.spread_volatility.abs() * s.power_forwards[i].abs() * t.sqrt())
            .max_val(T::from_f64(1e-6));
        let d_full = s_full / sig_t;
        let call_full = s_full * d_full.norm_cdf() + sig_t * d_full.norm_pdf();
        friction += call_full * mwh * df;
        derated += call_full * mwh * p_avail_t * df;
        let s_part =
            s.power_forwards[i] - s.gas_forwards[i] * s.part_load_heat_rate - s.vom_cost_per_mwh;
        let d_part = s_part / sig_t;
        let call_part = s_part * d_part.norm_cdf() + sig_t * d_part.norm_pdf();
        let res_price = s
            .spinning_reserve_prices
            .get(i)
            .copied()
            .unwrap_or_else(T::zero);
        let f_min = s
            .min_stable_load_fraction
            .max_val(T::zero())
            .min_val(T::one());
        let uplift_i = res_price * mwh * (T::one() - f_min) * p_avail_t * df;
        reserve += uplift_i;
        let pay_part = (f_min * call_part + (T::one() - f_min) * res_price) * mwh * p_avail_t * df;
        let pay_full = call_full * mwh * p_avail_t * df;
        if pay_part.val() >= pay_full.val() {
            opt += 1;
        }
        let c_start = s.startup_cost * s.capacity_mw * df;
        let c_ramp = s.fast_ramp_penalty_cost * s.capacity_mw * df;
        let next_0 = v_dp[0]
            .max_val(pay_part - c_start + v_dp[1])
            .max_val(pay_full - c_start - c_ramp + v_dp[2]);
        let next_1 = v_dp[0]
            .max_val(pay_part + v_dp[1])
            .max_val(pay_full - c_ramp + v_dp[2]);
        let next_2 = v_dp[0]
            .max_val(pay_part + v_dp[1])
            .max_val(pay_full + v_dp[2]);
        v_dp = [next_0, next_1, next_2];
    }
    let penalty = (friction - derated).max_val(T::zero());
    let constrained = v_dp[0].min_val(derated + reserve);
    EydelandWolyniecRampOutageVtaReport {
        steady_state_availability: avail,
        frictionless_spark_spread_npv: friction,
        reliability_derated_vta_npv: derated.min_val(friction),
        forced_outage_penalty_cost: penalty,
        ramp_and_startup_constrained_vta_npv: constrained,
        ancillary_spinning_reserve_uplift: reserve.max_val(T::zero()),
        optimal_part_load_reserve_periods: opt,
    }
}

#[cfg(test)]
mod domain_suite_35_tests {
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
mod domain_suite_26_tests {
    use super::*;
    use ql_ad::Dual;
    type D = Dual<1>;
    #[test]
    fn suite26_three_leg_weight3_ad() {
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
mod domain_suite_28_tests {
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
mod domain_suite_30_tests {
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
mod domain_suite_33_tests {
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

#[cfg(test)]
mod domain_suite_38_tests {
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

mod evaluate_eydeland_gibson_schwartz_and_kirk_spread_module {
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
            - (T::from_f64(2.0)
                * s.spot_cy_correlation
                * s.spot_volatility
                * s.convenience_yield_vol
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
    mod domain_suite_47_tests {
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
        fn suite47_commodities_invariants_and_dual_fd() {
            let r = evaluate_eydeland_gibson_schwartz_and_kirk_spread(&spec(0.0));
            assert!(r.kirk_spark_spread_call_pv >= r.intrinsic_spark_spread_pv - 1e-6);
            assert!(r.integrated_futures_variance > 0.0);
            let h = 1e-6;
            let fd = (evaluate_eydeland_gibson_schwartz_and_kirk_spread(&spec(h))
                .kirk_spark_spread_call_pv
                - evaluate_eydeland_gibson_schwartz_and_kirk_spread(&spec(-h))
                    .kirk_spark_spread_call_pv)
                / (2.0 * h);
            let ad = evaluate_eydeland_gibson_schwartz_and_kirk_spread(&spec(Dual::<1>::variable(
                0.0, 0,
            )))
            .kirk_spark_spread_call_pv;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_eydeland_gibson_schwartz_and_kirk_spread_module::*;

mod evaluate_eydeland_clean_spark_spread_option_module {
    // Eydeland & Wolyniec (2003) & Geman (2005) Power/Gas Clean Spark & Tolling Spread Option
    // via Bjerksund–Stensland (2014) / Margrabe Three-Leg Engine.

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
    fn price_spark<T: Real>(
        fp: T,
        fuel_leg: T,
        k_vom: T,
        sp: T,
        sg: T,
        rho: T,
        df: T,
        t: T,
    ) -> (T, T) {
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
    mod domain_suite_52_tests {
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
            let fd = (evaluate_eydeland_clean_spark_spread_option(&sp)
                .bjerksund_stensland_spark_call
                - evaluate_eydeland_clean_spark_spread_option(&sm).bjerksund_stensland_spark_call)
                / (2.0 * h);
            assert!((r_ad.bjerksund_stensland_spark_call.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_eydeland_clean_spark_spread_option_module::*;

mod evaluate_eydeland_crack_spread_three_asset_kirk_module {
    // Eydeland & Wolyniec (2003) & Carmona–Durrleman (2003) Three-Asset 3:2:1 Refinery
    // Crack-Spread Kirk–Margrabe Moment-Matching Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct EydelandThreeAssetCrackSpreadSpec<T: Real = f64> {
        pub crude_forward_barrel_fc: T,
        pub gasoline_forward_barrel_fg: T,
        pub distillate_forward_barrel_fh: T,
        pub fixed_refining_strike_k: T,
        pub crude_vol_sigma_c: T,
        pub gasoline_vol_sigma_g: T,
        pub distillate_vol_sigma_h: T,
        pub corr_gas_dist_rho_gh: T,
        pub corr_product_crude_rho_pc: T,
        pub risk_free_rate_r: T,
        pub expiry_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct EydelandThreeAssetCrackSpreadReport<T: Real = f64> {
        pub composite_product_basket_forward: T,
        pub intrinsic_321_crack_margin: T,
        pub composite_product_volatility: T,
        pub kirk_composite_crack_volatility: T,
        pub crack_spread_call_price_per_bbl: T,
        pub crude_leg_delta_sensitivity: T,
    }

    pub fn evaluate_eydeland_crack_spread_three_asset_kirk<T: Real>(
        spec: &EydelandThreeAssetCrackSpreadSpec<T>,
    ) -> EydelandThreeAssetCrackSpreadReport<T> {
        let fc = spec.crude_forward_barrel_fc.max_val(T::from_f64(1e-4));
        let fg = spec.gasoline_forward_barrel_fg.max_val(T::from_f64(1e-4));
        let fh = spec.distillate_forward_barrel_fh.max_val(T::from_f64(1e-4));
        let k = spec.fixed_refining_strike_k;
        let t = spec.expiry_years_t.max_val(T::from_f64(0.02));
        let df = (-spec.risk_free_rate_r * t).exp();

        // 3:2:1 crack per barrel of crude: ProductBasket = (2/3) * F_gas + (1/3) * F_distillate:
        let w_g = T::from_f64(2.0 / 3.0);
        let w_h = T::from_f64(1.0 / 3.0);
        let f_prod = w_g * fg + w_h * fh;
        let intrinsic = f_prod - fc - k;

        // Levy (1992) moment-matched volatility of the 2:1 product basket:
        let share_g = (w_g * fg) / f_prod;
        let share_h = (w_h * fh) / f_prod;
        let var_prod = (share_g * share_g * spec.gasoline_vol_sigma_g * spec.gasoline_vol_sigma_g
            + share_h * share_h * spec.distillate_vol_sigma_h * spec.distillate_vol_sigma_h
            + T::from_f64(2.0)
                * share_g
                * share_h
                * spec.corr_gas_dist_rho_gh
                * spec.gasoline_vol_sigma_g
                * spec.distillate_vol_sigma_h)
            .max_val(T::from_f64(1e-8));
        let sig_prod = var_prod.sqrt();

        // Kirk (1995) approximation against composite strike A = F_crude + K:
        let a = (fc + k).max_val(T::from_f64(1e-4));
        let w_c = fc / a;
        let var_kirk = (var_prod
            - T::from_f64(2.0)
                * spec.corr_product_crude_rho_pc
                * sig_prod
                * spec.crude_vol_sigma_c
                * w_c
            + spec.crude_vol_sigma_c * spec.crude_vol_sigma_c * w_c * w_c)
            .max_val(T::from_f64(1e-8));
        let sig_kirk = var_kirk.sqrt();
        let std = sig_kirk * t.sqrt();

        let d1 = ((f_prod / a).ln() + T::from_f64(0.5) * std * std) / std;
        let d2 = d1 - std;
        let call_price = df * (f_prod * norm_cdf(d1) - a * norm_cdf(d2));
        let crude_delta = -df * norm_cdf(d2);

        EydelandThreeAssetCrackSpreadReport {
            composite_product_basket_forward: f_prod,
            intrinsic_321_crack_margin: intrinsic,
            composite_product_volatility: sig_prod,
            kirk_composite_crack_volatility: sig_kirk,
            crack_spread_call_price_per_bbl: call_price,
            crude_leg_delta_sensitivity: crude_delta,
        }
    }

    #[cfg(test)]
    mod domain_suite_53_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn eydeland_crack_spread_three_asset_kirk_invariants_and_ad_parity() {
            let spec = EydelandThreeAssetCrackSpreadSpec {
                crude_forward_barrel_fc: 76.0,
                gasoline_forward_barrel_fg: 94.0,
                distillate_forward_barrel_fh: 98.0,
                fixed_refining_strike_k: 6.5,
                crude_vol_sigma_c: 0.34,
                gasoline_vol_sigma_g: 0.38,
                distillate_vol_sigma_h: 0.36,
                corr_gas_dist_rho_gh: 0.82,
                corr_product_crude_rho_pc: 0.86,
                risk_free_rate_r: 0.04,
                expiry_years_t: 0.5,
            };
            let rep = evaluate_eydeland_crack_spread_three_asset_kirk(&spec);
            assert!(rep.intrinsic_321_crack_margin > 0.0);
            assert!(rep.crack_spread_call_price_per_bbl > rep.intrinsic_321_crack_margin * 0.9);
            assert!(rep.crude_leg_delta_sensitivity < 0.0);

            let h = 1e-6;
            let r_ad = evaluate_eydeland_crack_spread_three_asset_kirk(
                &EydelandThreeAssetCrackSpreadSpec {
                    crude_forward_barrel_fc: Dual::<1>::variable(spec.crude_forward_barrel_fc, 0),
                    gasoline_forward_barrel_fg: Dual::<1>::from_value(
                        spec.gasoline_forward_barrel_fg,
                    ),
                    distillate_forward_barrel_fh: Dual::<1>::from_value(
                        spec.distillate_forward_barrel_fh,
                    ),
                    fixed_refining_strike_k: Dual::<1>::from_value(spec.fixed_refining_strike_k),
                    crude_vol_sigma_c: Dual::<1>::from_value(spec.crude_vol_sigma_c),
                    gasoline_vol_sigma_g: Dual::<1>::from_value(spec.gasoline_vol_sigma_g),
                    distillate_vol_sigma_h: Dual::<1>::from_value(spec.distillate_vol_sigma_h),
                    corr_gas_dist_rho_gh: Dual::<1>::from_value(spec.corr_gas_dist_rho_gh),
                    corr_product_crude_rho_pc: Dual::<1>::from_value(
                        spec.corr_product_crude_rho_pc,
                    ),
                    risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
                    expiry_years_t: Dual::<1>::from_value(spec.expiry_years_t),
                },
            );
            let mut sp = spec;
            sp.crude_forward_barrel_fc += h;
            let mut sm = spec;
            sm.crude_forward_barrel_fc -= h;
            let fd = (evaluate_eydeland_crack_spread_three_asset_kirk(&sp)
                .crack_spread_call_price_per_bbl
                - evaluate_eydeland_crack_spread_three_asset_kirk(&sm)
                    .crack_spread_call_price_per_bbl)
                / (2.0 * h);
            assert!((r_ad.crack_spread_call_price_per_bbl.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_eydeland_crack_spread_three_asset_kirk_module::*;

mod evaluate_bessembinder_lemmon_power_forward_premium_module {
    // Geman (2005) & Bessembinder–Lemmon (2002) Equilibrium Electricity Forward Risk
    // Premium (`Variance vs. Positive Spike Skewness`) & Reserve-Margin Convexity Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BessembinderLemmonPowerForwardPremiumSpec<T: Real = f64> {
        pub baseline_marginal_fuel_cost_c0: T,
        pub expected_load_to_capacity_ratio_u: T,
        pub load_shock_volatility_sigma_u: T,
        pub supply_stack_convexity_exponent_gamma: T,
        pub peaking_merit_slope_alpha: T,
        pub producer_variance_aversion_a_var: T,
        pub retailer_spike_skew_preference_b_skew: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BessembinderLemmonPowerForwardPremiumReport<T: Real = f64> {
        pub expected_spot_electricity_price: T,
        pub spot_price_variance: T,
        pub spot_price_cubic_skewness_moment: T,
        pub variance_hedging_discount_component: T,
        pub spike_skewness_premium_component: T,
        pub net_equilibrium_forward_risk_premium: T,
        pub equilibrium_power_forward_price: T,
        pub reserve_margin_tightness_index: T,
    }

    pub fn evaluate_bessembinder_lemmon_power_forward_premium<T: Real>(
        spec: &BessembinderLemmonPowerForwardPremiumSpec<T>,
    ) -> BessembinderLemmonPowerForwardPremiumReport<T> {
        let c0 = spec
            .baseline_marginal_fuel_cost_c0
            .max_val(T::from_f64(5.0));
        let u = spec
            .expected_load_to_capacity_ratio_u
            .max_val(T::from_f64(0.20))
            .min_val(T::from_f64(0.98));
        let sig_u = spec
            .load_shock_volatility_sigma_u
            .max_val(T::from_f64(0.01));
        let gam = spec
            .supply_stack_convexity_exponent_gamma
            .max_val(T::from_f64(1.5));
        let alpha = spec.peaking_merit_slope_alpha.max_val(T::from_f64(1.0));
        let a_var = spec.producer_variance_aversion_a_var.max_val(T::zero());
        let b_skew = spec
            .retailer_spike_skew_preference_b_skew
            .max_val(T::zero());

        // Convex generation stack P(U) = c0 + alpha * U^gamma:
        // Using Taylor-Itô expansion around expected utilization u:
        let u_pow_g = (gam * u.ln()).exp();
        let p_at_mean = c0 + alpha * u_pow_g;
        let p1 = alpha * gam * (u_pow_g / u);
        let p2 = alpha * gam * (gam - T::one()) * (u_pow_g / (u * u));

        let exp_spot = p_at_mean + T::from_f64(0.5) * p2 * sig_u * sig_u;
        let var_spot =
            p1 * p1 * sig_u * sig_u + T::from_f64(0.5) * p2 * p2 * sig_u * sig_u * sig_u * sig_u;
        // Third central moment E[(S - E[S])^3] = 3 p1^2 p2 sigma_u^4 + p2^3 sigma_u^6 > 0 for convex stack (gamma > 1):
        let skew3_spot = T::from_f64(3.0) * p1 * p1 * p2 * sig_u * sig_u * sig_u * sig_u
            + p2 * p2 * p2 * sig_u * sig_u * sig_u * sig_u * sig_u * sig_u;

        // Bessembinder & Lemmon (2002, JF Eq. 8): F - E[S] = -a_var * Var(S) + b_skew * Skew_3(S)
        let var_comp = -a_var * var_spot;
        let skew_comp = b_skew * skew3_spot;
        let net_prem = var_comp + skew_comp;
        let fwd_price = exp_spot + net_prem;
        let tightness = u / (T::one() - u).max_val(T::from_f64(0.02));

        BessembinderLemmonPowerForwardPremiumReport {
            expected_spot_electricity_price: exp_spot,
            spot_price_variance: var_spot,
            spot_price_cubic_skewness_moment: skew3_spot,
            variance_hedging_discount_component: var_comp,
            spike_skewness_premium_component: skew_comp,
            net_equilibrium_forward_risk_premium: net_prem,
            equilibrium_power_forward_price: fwd_price,
            reserve_margin_tightness_index: tightness,
        }
    }

    #[cfg(test)]
    mod evaluate_bessembinder_lemmon_power_forward_premium_module_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn bessembinder_lemmon_power_forward_premium_invariants_and_ad_parity() {
            let spec = BessembinderLemmonPowerForwardPremiumSpec {
                baseline_marginal_fuel_cost_c0: 28.0,
                expected_load_to_capacity_ratio_u: 0.84,
                load_shock_volatility_sigma_u: 0.12,
                supply_stack_convexity_exponent_gamma: 4.5,
                peaking_merit_slope_alpha: 65.0,
                producer_variance_aversion_a_var: 0.015,
                retailer_spike_skew_preference_b_skew: 0.0045,
            };
            let rep = evaluate_bessembinder_lemmon_power_forward_premium(&spec);
            assert!(rep.spot_price_cubic_skewness_moment > 0.0);
            assert!(rep.net_equilibrium_forward_risk_premium > 0.0);
            assert!(rep.equilibrium_power_forward_price > rep.expected_spot_electricity_price);

            let h = 1e-6;
            let r_ad = evaluate_bessembinder_lemmon_power_forward_premium(
                &BessembinderLemmonPowerForwardPremiumSpec {
                    baseline_marginal_fuel_cost_c0: Dual::<7>::variable(
                        spec.baseline_marginal_fuel_cost_c0,
                        0,
                    ),
                    expected_load_to_capacity_ratio_u: Dual::<7>::variable(
                        spec.expected_load_to_capacity_ratio_u,
                        1,
                    ),
                    load_shock_volatility_sigma_u: Dual::<7>::variable(
                        spec.load_shock_volatility_sigma_u,
                        2,
                    ),
                    supply_stack_convexity_exponent_gamma: Dual::<7>::variable(
                        spec.supply_stack_convexity_exponent_gamma,
                        3,
                    ),
                    peaking_merit_slope_alpha: Dual::<7>::variable(
                        spec.peaking_merit_slope_alpha,
                        4,
                    ),
                    producer_variance_aversion_a_var: Dual::<7>::variable(
                        spec.producer_variance_aversion_a_var,
                        5,
                    ),
                    retailer_spike_skew_preference_b_skew: Dual::<7>::variable(
                        spec.retailer_spike_skew_preference_b_skew,
                        6,
                    ),
                },
            );
            for (idx, &ad_deriv) in r_ad.equilibrium_power_forward_price.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.baseline_marginal_fuel_cost_c0 += h;
                        sm.baseline_marginal_fuel_cost_c0 -= h;
                    }
                    1 => {
                        sp.expected_load_to_capacity_ratio_u += h;
                        sm.expected_load_to_capacity_ratio_u -= h;
                    }
                    2 => {
                        sp.load_shock_volatility_sigma_u += h;
                        sm.load_shock_volatility_sigma_u -= h;
                    }
                    3 => {
                        sp.supply_stack_convexity_exponent_gamma += h;
                        sm.supply_stack_convexity_exponent_gamma -= h;
                    }
                    4 => {
                        sp.peaking_merit_slope_alpha += h;
                        sm.peaking_merit_slope_alpha -= h;
                    }
                    5 => {
                        sp.producer_variance_aversion_a_var += h;
                        sm.producer_variance_aversion_a_var -= h;
                    }
                    _ => {
                        sp.retailer_spike_skew_preference_b_skew += h;
                        sm.retailer_spike_skew_preference_b_skew -= h;
                    }
                }
                let fd = (evaluate_bessembinder_lemmon_power_forward_premium(&sp)
                    .equilibrium_power_forward_price
                    - evaluate_bessembinder_lemmon_power_forward_premium(&sm)
                        .equilibrium_power_forward_price)
                    / (2.0 * h);
                assert!((ad_deriv - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_bessembinder_lemmon_power_forward_premium_module::*;

mod evaluate_deng_johnson_locational_transmission_congestion_option_module {
    // Deng, Johnson & Sogomonian (2001, Energy Economics) & Eydeland-Wolyniec (2003, Ch. 9)
    // Locational Electricity Transmission Congestion Option (FTR), Line-Loss Factor & Congestion Rent Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct DengJohnsonTransmissionCongestionOptionSpec<T: Real = f64> {
        pub sink_hub_forward_price_fb: T,
        pub source_hub_forward_price_fa: T,
        pub transmission_loss_factor_eta: T,
        pub wheeling_tariff_strike_k_trans: T,
        pub sink_hub_volatility_sigma_b: T,
        pub source_hub_volatility_sigma_a: T,
        pub hub_price_correlation_rho_ab: T,
        pub risk_free_rate_r: T,
        pub expiry_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct DengJohnsonTransmissionCongestionOptionReport<T: Real = f64> {
        pub loss_adjusted_sink_forward: T,
        pub kirk_effective_spread_volatility: T,
        pub ftr_option_value_per_mwh: T,
        pub ftr_obligation_value_per_mwh: T,
        pub counterflow_optionality_value_per_mwh: T,
        pub sink_hub_shadow_delta: T,
    }

    pub fn evaluate_deng_johnson_locational_transmission_congestion_option<T: Real>(
        spec: &DengJohnsonTransmissionCongestionOptionSpec<T>,
    ) -> DengJohnsonTransmissionCongestionOptionReport<T> {
        let fb = spec.sink_hub_forward_price_fb.max_val(T::from_f64(1.0));
        let fa = spec.source_hub_forward_price_fa.max_val(T::from_f64(1.0));
        let eta = spec
            .transmission_loss_factor_eta
            .max_val(T::from_f64(0.50))
            .min_val(T::one());
        let k_tr = spec
            .wheeling_tariff_strike_k_trans
            .max_val(T::from_f64(0.1));
        let sig_b = spec.sink_hub_volatility_sigma_b.max_val(T::from_f64(0.05));
        let sig_a = spec
            .source_hub_volatility_sigma_a
            .max_val(T::from_f64(0.05));
        let rho = spec
            .hub_price_correlation_rho_ab
            .max_val(T::from_f64(-0.95))
            .min_val(T::from_f64(0.95));
        let r = spec.risk_free_rate_r;
        let t = spec.expiry_years_t.max_val(T::from_f64(0.05));

        let fb_eff = eta * fb;
        let f_comp = fa + k_tr;
        let w_a = fa / f_comp;

        let var_kirk = (sig_b * sig_b - T::from_f64(2.0) * w_a * rho * sig_a * sig_b
            + w_a * w_a * sig_a * sig_a)
            .max_val(T::from_f64(1e-4));
        let sig_kirk = var_kirk.sqrt();
        let std = sig_kirk * t.sqrt();

        let d1 = ((fb_eff / f_comp).ln() + T::from_f64(0.5) * var_kirk * t) / std;
        let d2 = d1 - std;

        let df = (-r * t).exp();
        let ftr_opt = df * (fb_eff * norm_cdf(d1) - f_comp * norm_cdf(d2));
        let ftr_ob = df * (fb_eff - f_comp);
        let counterflow_opt = ftr_opt - ftr_ob;
        let sink_delta = eta * df * norm_cdf(d1);

        DengJohnsonTransmissionCongestionOptionReport {
            loss_adjusted_sink_forward: fb_eff,
            kirk_effective_spread_volatility: sig_kirk,
            ftr_option_value_per_mwh: ftr_opt,
            ftr_obligation_value_per_mwh: ftr_ob,
            counterflow_optionality_value_per_mwh: counterflow_opt,
            sink_hub_shadow_delta: sink_delta,
        }
    }

    #[cfg(test)]
    mod evaluate_deng_johnson_locational_transmission_congestion_option_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn deng_johnson_locational_transmission_option_invariants_and_ad_parity() {
            let spec = DengJohnsonTransmissionCongestionOptionSpec {
                sink_hub_forward_price_fb: 68.0,
                source_hub_forward_price_fa: 54.0,
                transmission_loss_factor_eta: 0.96,
                wheeling_tariff_strike_k_trans: 4.5,
                sink_hub_volatility_sigma_b: 0.48,
                source_hub_volatility_sigma_a: 0.42,
                hub_price_correlation_rho_ab: 0.78,
                risk_free_rate_r: 0.04,
                expiry_years_t: 0.50,
            };
            let rep = evaluate_deng_johnson_locational_transmission_congestion_option(&spec);
            assert!(rep.ftr_option_value_per_mwh > rep.ftr_obligation_value_per_mwh);
            assert!(rep.counterflow_optionality_value_per_mwh > 0.0);
            assert!(rep.sink_hub_shadow_delta > 0.0);

            let h = 1e-5;
            let r_ad = evaluate_deng_johnson_locational_transmission_congestion_option(
                &DengJohnsonTransmissionCongestionOptionSpec {
                    sink_hub_forward_price_fb: Dual::<9>::variable(
                        spec.sink_hub_forward_price_fb,
                        0,
                    ),
                    source_hub_forward_price_fa: Dual::<9>::variable(
                        spec.source_hub_forward_price_fa,
                        1,
                    ),
                    transmission_loss_factor_eta: Dual::<9>::variable(
                        spec.transmission_loss_factor_eta,
                        2,
                    ),
                    wheeling_tariff_strike_k_trans: Dual::<9>::variable(
                        spec.wheeling_tariff_strike_k_trans,
                        3,
                    ),
                    sink_hub_volatility_sigma_b: Dual::<9>::variable(
                        spec.sink_hub_volatility_sigma_b,
                        4,
                    ),
                    source_hub_volatility_sigma_a: Dual::<9>::variable(
                        spec.source_hub_volatility_sigma_a,
                        5,
                    ),
                    hub_price_correlation_rho_ab: Dual::<9>::variable(
                        spec.hub_price_correlation_rho_ab,
                        6,
                    ),
                    risk_free_rate_r: Dual::<9>::variable(spec.risk_free_rate_r, 7),
                    expiry_years_t: Dual::<9>::variable(spec.expiry_years_t, 8),
                },
            );
            for (idx, &ad) in r_ad.ftr_option_value_per_mwh.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.sink_hub_forward_price_fb += h;
                        sm.sink_hub_forward_price_fb -= h;
                    }
                    1 => {
                        sp.source_hub_forward_price_fa += h;
                        sm.source_hub_forward_price_fa -= h;
                    }
                    2 => {
                        sp.transmission_loss_factor_eta += h;
                        sm.transmission_loss_factor_eta -= h;
                    }
                    3 => {
                        sp.wheeling_tariff_strike_k_trans += h;
                        sm.wheeling_tariff_strike_k_trans -= h;
                    }
                    4 => {
                        sp.sink_hub_volatility_sigma_b += h;
                        sm.sink_hub_volatility_sigma_b -= h;
                    }
                    5 => {
                        sp.source_hub_volatility_sigma_a += h;
                        sm.source_hub_volatility_sigma_a -= h;
                    }
                    6 => {
                        sp.hub_price_correlation_rho_ab += h;
                        sm.hub_price_correlation_rho_ab -= h;
                    }
                    7 => {
                        sp.risk_free_rate_r += h;
                        sm.risk_free_rate_r -= h;
                    }
                    _ => {
                        sp.expiry_years_t += h;
                        sm.expiry_years_t -= h;
                    }
                }
                let fd = (evaluate_deng_johnson_locational_transmission_congestion_option(&sp)
                    .ftr_option_value_per_mwh
                    - evaluate_deng_johnson_locational_transmission_congestion_option(&sm)
                        .ftr_option_value_per_mwh)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_deng_johnson_locational_transmission_congestion_option_module::*;

mod evaluate_carmona_durrleman_spark_spread_lower_bound_module {
    // Carmona & Durrleman (2003, SIAM Review 45(4), "Pricing and Hedging Spread Options"),
    // Deng, Johnson & Sogomonian (2001, Review of Derivatives Research 4(3)),
    // & Eydeland & Wolyniec (2003, Energy and Power Risk Management, Ch. 7)
    // Power Plant Tolling Spark Spread Lower-Bound & Second-Order Curvature Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct CarmonaDurrlemanSparkSpreadSpec<T: Real = f64> {
        pub electricity_forward_f_e: T,
        pub natural_gas_forward_f_g: T,
        pub heat_rate_mmbtu_per_mwh_h: T,
        pub variable_om_strike_k_om: T,
        pub electricity_volatility_sigma_e: T,
        pub gas_volatility_sigma_g: T,
        pub power_gas_correlation_rho: T,
        pub risk_free_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct CarmonaDurrlemanSparkSpreadReport<T: Real = f64> {
        pub kirk_composite_volatility: T,
        pub carmona_durrleman_effective_volatility: T,
        pub kirk_benchmark_spark_call_value: T,
        pub carmona_durrleman_spark_call_value: T,
        pub curvature_basis_over_kirk: T,
        pub power_forward_delta: T,
        pub implied_plant_dispatch_probability: T,
    }

    pub fn evaluate_carmona_durrleman_spark_spread_lower_bound<T: Real>(
        spec: &CarmonaDurrlemanSparkSpreadSpec<T>,
    ) -> CarmonaDurrlemanSparkSpreadReport<T> {
        let f_e = spec.electricity_forward_f_e.max_val(T::from_f64(1.0));
        let f_g = spec.natural_gas_forward_f_g.max_val(T::from_f64(0.1));
        let h = spec.heat_rate_mmbtu_per_mwh_h.max_val(T::from_f64(1.0));
        let k_om = spec.variable_om_strike_k_om.max_val(T::from_f64(0.1));
        let sig_e = spec
            .electricity_volatility_sigma_e
            .max_val(T::from_f64(0.02));
        let sig_g = spec.gas_volatility_sigma_g.max_val(T::from_f64(0.02));
        let rho = spec
            .power_gas_correlation_rho
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));
        let r = spec.risk_free_rate_r;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let x_fuel = h * f_g;
        let x_tot = x_fuel + k_om;
        let psi = x_fuel / x_tot;
        let om_share = T::one() - psi;

        let var_kirk = (sig_e * sig_e - T::from_f64(2.0) * rho * psi * sig_e * sig_g
            + psi * psi * sig_g * sig_g)
            .max_val(T::from_f64(1e-6));
        let sig_kirk = var_kirk.sqrt();

        // Carmona-Durrleman (2003) / Bjerksund-Stensland second-order log-sum curvature variance correction:
        let curv_var = T::from_f64(0.25)
            * (T::one() - rho * rho)
            * psi
            * psi
            * om_share
            * om_share
            * (sig_g * sig_g * sig_g * sig_g / var_kirk)
            * t;
        let sig_cd = (var_kirk + curv_var).sqrt();

        let sqrt_t = t.sqrt();
        let df = (-r * t).exp();
        let log_m = (f_e / x_tot).ln();

        let std_kirk = sig_kirk * sqrt_t;
        let d1_k = (log_m + T::from_f64(0.5) * var_kirk * t) / std_kirk;
        let d2_k = d1_k - std_kirk;
        let c_kirk = df * (f_e * norm_cdf(d1_k) - x_tot * norm_cdf(d2_k));

        let std_cd = sig_cd * sqrt_t;
        let d1_cd = (log_m + T::from_f64(0.5) * sig_cd * sig_cd * t) / std_cd;
        let d2_cd = d1_cd - std_cd;
        let c_cd = df * (f_e * norm_cdf(d1_cd) - x_tot * norm_cdf(d2_cd));

        let basis = c_cd - c_kirk;
        let delta_e = df * norm_cdf(d1_cd);
        let dispatch_prob = norm_cdf(d2_cd);

        CarmonaDurrlemanSparkSpreadReport {
            kirk_composite_volatility: sig_kirk,
            carmona_durrleman_effective_volatility: sig_cd,
            kirk_benchmark_spark_call_value: c_kirk,
            carmona_durrleman_spark_call_value: c_cd,
            curvature_basis_over_kirk: basis,
            power_forward_delta: delta_e,
            implied_plant_dispatch_probability: dispatch_prob,
        }
    }

    #[cfg(test)]
    mod evaluate_carmona_durrleman_spark_spread_lower_bound_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn carmona_durrleman_spark_spread_invariants_and_ad_parity() {
            let spec = CarmonaDurrlemanSparkSpreadSpec {
                electricity_forward_f_e: 62.0,
                natural_gas_forward_f_g: 7.20,
                heat_rate_mmbtu_per_mwh_h: 7.50,
                variable_om_strike_k_om: 6.50,
                electricity_volatility_sigma_e: 0.48,
                gas_volatility_sigma_g: 0.42,
                power_gas_correlation_rho: 0.65,
                risk_free_rate_r: 0.04,
                maturity_years_t: 0.75,
            };
            let rep = evaluate_carmona_durrleman_spark_spread_lower_bound(&spec);
            assert!(rep.carmona_durrleman_effective_volatility > rep.kirk_composite_volatility);
            assert!(rep.curvature_basis_over_kirk > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_carmona_durrleman_spark_spread_lower_bound(
                &CarmonaDurrlemanSparkSpreadSpec {
                    electricity_forward_f_e: Dual::<6>::variable(spec.electricity_forward_f_e, 0),
                    natural_gas_forward_f_g: Dual::<6>::variable(spec.natural_gas_forward_f_g, 1),
                    heat_rate_mmbtu_per_mwh_h: Dual::from_f64(spec.heat_rate_mmbtu_per_mwh_h),
                    variable_om_strike_k_om: Dual::<6>::variable(spec.variable_om_strike_k_om, 2),
                    electricity_volatility_sigma_e: Dual::<6>::variable(
                        spec.electricity_volatility_sigma_e,
                        3,
                    ),
                    gas_volatility_sigma_g: Dual::<6>::variable(spec.gas_volatility_sigma_g, 4),
                    power_gas_correlation_rho: Dual::<6>::variable(
                        spec.power_gas_correlation_rho,
                        5,
                    ),
                    risk_free_rate_r: Dual::from_f64(spec.risk_free_rate_r),
                    maturity_years_t: Dual::from_f64(spec.maturity_years_t),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.electricity_forward_f_e,
                    spec.natural_gas_forward_f_g,
                    spec.variable_om_strike_k_om,
                    spec.electricity_volatility_sigma_e,
                    spec.gas_volatility_sigma_g,
                    spec.power_gas_correlation_rho,
                ];
                v[i] += h;
                let up = evaluate_carmona_durrleman_spark_spread_lower_bound(
                    &CarmonaDurrlemanSparkSpreadSpec {
                        electricity_forward_f_e: v[0],
                        natural_gas_forward_f_g: v[1],
                        variable_om_strike_k_om: v[2],
                        electricity_volatility_sigma_e: v[3],
                        gas_volatility_sigma_g: v[4],
                        power_gas_correlation_rho: v[5],
                        ..spec
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_carmona_durrleman_spark_spread_lower_bound(
                    &CarmonaDurrlemanSparkSpreadSpec {
                        electricity_forward_f_e: v[0],
                        natural_gas_forward_f_g: v[1],
                        variable_om_strike_k_om: v[2],
                        electricity_volatility_sigma_e: v[3],
                        gas_volatility_sigma_g: v[4],
                        power_gas_correlation_rho: v[5],
                        ..spec
                    },
                );
                let fd = (up.carmona_durrleman_spark_call_value
                    - dn.carmona_durrleman_spark_call_value)
                    / (2.0 * h);
                assert!(
                    (r_ad.carmona_durrleman_spark_call_value.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_carmona_durrleman_spark_spread_lower_bound_module::*;

mod evaluate_eydeland_wolyniec_virtual_power_plant_heat_rate_strip_module {
    // Eydeland & Wolyniec (2003, Energy and Power Risk Management, Ch. 8,
    // "Power Plant Valuation and Tolling Agreements"), Kirk (1995), & Carmona & Durrleman (2003)
    // Virtual Power Plant (VPP) Heat-Rate Tolling Strip & Startup Friction Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct EydelandWolyniecVppHeatRateStripSpec<T: Real = f64> {
        pub power_forward_price_f_p: T,
        pub gas_forward_price_f_g: T,
        pub plant_heat_rate_h: T,
        pub variable_om_cost_v_om: T,
        pub startup_cost_per_mwh_c_start: T,
        pub power_volatility_sigma_p: T,
        pub gas_volatility_sigma_g: T,
        pub power_gas_correlation_rho: T,
        pub risk_free_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct EydelandWolyniecVppHeatRateStripReport<T: Real = f64> {
        pub effective_all_in_dispatch_strike: T,
        pub kirk_composite_spread_volatility: T,
        pub frictionless_spark_spread_call_value: T,
        pub vpp_startup_adjusted_tolling_value: T,
        pub startup_cost_option_value_drag: T,
        pub expected_dispatch_capacity_factor: T,
    }

    pub fn evaluate_eydeland_wolyniec_virtual_power_plant_heat_rate_strip<T: Real>(
        spec: &EydelandWolyniecVppHeatRateStripSpec<T>,
    ) -> EydelandWolyniecVppHeatRateStripReport<T> {
        let f_p = spec.power_forward_price_f_p.max_val(T::from_f64(5.0));
        let f_g = spec.gas_forward_price_f_g.max_val(T::from_f64(0.5));
        let h_rate = spec.plant_heat_rate_h.max_val(T::from_f64(4.0));
        let v_om = spec.variable_om_cost_v_om.max_val(T::from_f64(0.1));
        let c_start = spec.startup_cost_per_mwh_c_start.max_val(T::from_f64(0.01));
        let sig_p = spec.power_volatility_sigma_p.max_val(T::from_f64(0.02));
        let sig_g = spec.gas_volatility_sigma_g.max_val(T::from_f64(0.02));
        let rho = spec
            .power_gas_correlation_rho
            .max_val(T::from_f64(-0.95))
            .min_val(T::from_f64(0.95));
        let r = spec.risk_free_rate_r.max_val(T::from_f64(0.001));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let fuel_cost = h_rate * f_g;
        let k_spark = fuel_cost + v_om;
        let k_vpp = k_spark + c_start;

        let df = (-r * t).exp();
        let sqrt_t = t.sqrt();

        // Kirk (1995) / Eydeland-Wolyniec (2003, Ch. 8) composite volatility for frictionless spark spread:
        let w_g_spark = (fuel_cost / k_spark) * sig_g;
        let sig_kirk_spark = (sig_p * sig_p - T::from_f64(2.0) * rho * sig_p * w_g_spark
            + w_g_spark * w_g_spark)
            .max_val(T::from_f64(1e-4))
            .sqrt();
        let std_spark = sig_kirk_spark * sqrt_t;
        let d1_spark = ((f_p / k_spark).ln()
            + T::from_f64(0.5) * sig_kirk_spark * sig_kirk_spark * t)
            / std_spark;
        let d2_spark = d1_spark - std_spark;
        let val_spark = df * (f_p * norm_cdf(d1_spark) - k_spark * norm_cdf(d2_spark));

        // Kirk composite volatility for VPP dispatch block including amortized startup cost C_start:
        let w_g_vpp = (fuel_cost / k_vpp) * sig_g;
        let sig_kirk_vpp = (sig_p * sig_p - T::from_f64(2.0) * rho * sig_p * w_g_vpp
            + w_g_vpp * w_g_vpp)
            .max_val(T::from_f64(1e-4))
            .sqrt();
        let std_vpp = sig_kirk_vpp * sqrt_t;
        let d1_vpp =
            ((f_p / k_vpp).ln() + T::from_f64(0.5) * sig_kirk_vpp * sig_kirk_vpp * t) / std_vpp;
        let d2_vpp = d1_vpp - std_vpp;
        let val_vpp = df * (f_p * norm_cdf(d1_vpp) - k_vpp * norm_cdf(d2_vpp));

        let drag = val_spark - val_vpp;
        let cap_factor = norm_cdf(d2_vpp);

        EydelandWolyniecVppHeatRateStripReport {
            effective_all_in_dispatch_strike: k_vpp,
            kirk_composite_spread_volatility: sig_kirk_vpp,
            frictionless_spark_spread_call_value: val_spark,
            vpp_startup_adjusted_tolling_value: val_vpp,
            startup_cost_option_value_drag: drag,
            expected_dispatch_capacity_factor: cap_factor,
        }
    }

    #[cfg(test)]
    mod evaluate_eydeland_wolyniec_virtual_power_plant_heat_rate_strip_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn eydeland_wolyniec_vpp_heat_rate_strip_invariants_and_ad_parity() {
            let spec = EydelandWolyniecVppHeatRateStripSpec {
                power_forward_price_f_p: 68.0,
                gas_forward_price_f_g: 7.50,
                plant_heat_rate_h: 7.80,
                variable_om_cost_v_om: 3.50,
                startup_cost_per_mwh_c_start: 2.20,
                power_volatility_sigma_p: 0.48,
                gas_volatility_sigma_g: 0.36,
                power_gas_correlation_rho: 0.72,
                risk_free_rate_r: 0.04,
                maturity_years_t: 0.75,
            };
            let rep = evaluate_eydeland_wolyniec_virtual_power_plant_heat_rate_strip(&spec);
            assert!(
                rep.frictionless_spark_spread_call_value > rep.vpp_startup_adjusted_tolling_value
            );
            assert!(rep.startup_cost_option_value_drag > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_eydeland_wolyniec_virtual_power_plant_heat_rate_strip(
                &EydelandWolyniecVppHeatRateStripSpec {
                    power_forward_price_f_p: Dual::<10>::variable(spec.power_forward_price_f_p, 0),
                    gas_forward_price_f_g: Dual::<10>::variable(spec.gas_forward_price_f_g, 1),
                    plant_heat_rate_h: Dual::<10>::variable(spec.plant_heat_rate_h, 2),
                    variable_om_cost_v_om: Dual::<10>::variable(spec.variable_om_cost_v_om, 3),
                    startup_cost_per_mwh_c_start: Dual::<10>::variable(
                        spec.startup_cost_per_mwh_c_start,
                        4,
                    ),
                    power_volatility_sigma_p: Dual::<10>::variable(
                        spec.power_volatility_sigma_p,
                        5,
                    ),
                    gas_volatility_sigma_g: Dual::<10>::variable(spec.gas_volatility_sigma_g, 6),
                    power_gas_correlation_rho: Dual::<10>::variable(
                        spec.power_gas_correlation_rho,
                        7,
                    ),
                    risk_free_rate_r: Dual::<10>::variable(spec.risk_free_rate_r, 8),
                    maturity_years_t: Dual::<10>::variable(spec.maturity_years_t, 9),
                },
            );
            for i in 0..10 {
                let mut v = [
                    spec.power_forward_price_f_p,
                    spec.gas_forward_price_f_g,
                    spec.plant_heat_rate_h,
                    spec.variable_om_cost_v_om,
                    spec.startup_cost_per_mwh_c_start,
                    spec.power_volatility_sigma_p,
                    spec.gas_volatility_sigma_g,
                    spec.power_gas_correlation_rho,
                    spec.risk_free_rate_r,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_eydeland_wolyniec_virtual_power_plant_heat_rate_strip(
                    &EydelandWolyniecVppHeatRateStripSpec {
                        power_forward_price_f_p: v[0],
                        gas_forward_price_f_g: v[1],
                        plant_heat_rate_h: v[2],
                        variable_om_cost_v_om: v[3],
                        startup_cost_per_mwh_c_start: v[4],
                        power_volatility_sigma_p: v[5],
                        gas_volatility_sigma_g: v[6],
                        power_gas_correlation_rho: v[7],
                        risk_free_rate_r: v[8],
                        maturity_years_t: v[9],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_eydeland_wolyniec_virtual_power_plant_heat_rate_strip(
                    &EydelandWolyniecVppHeatRateStripSpec {
                        power_forward_price_f_p: v[0],
                        gas_forward_price_f_g: v[1],
                        plant_heat_rate_h: v[2],
                        variable_om_cost_v_om: v[3],
                        startup_cost_per_mwh_c_start: v[4],
                        power_volatility_sigma_p: v[5],
                        gas_volatility_sigma_g: v[6],
                        power_gas_correlation_rho: v[7],
                        risk_free_rate_r: v[8],
                        maturity_years_t: v[9],
                    },
                );
                let fd = (up.vpp_startup_adjusted_tolling_value
                    - dn.vpp_startup_adjusted_tolling_value)
                    / (2.0 * h);
                assert!(
                    (r_ad.vpp_startup_adjusted_tolling_value.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_eydeland_wolyniec_virtual_power_plant_heat_rate_strip_module::*;
