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

#[derive(Copy,Clone,Debug,PartialEq,Eq)]
pub enum DegreeDayKind{Hdd,Cdd}
#[derive(Clone,Debug)]
pub struct WeatherDerivativeSpec<T: Real=f64>{pub kind:DegreeDayKind,pub temperature_base:T,pub initial_temperature:T,pub mean_level:T,pub seasonal_amplitude:T,pub seasonal_phase:T,pub mean_reversion:T,pub daily_vol:T,pub market_price_risk:T,pub times:Vec<f64>,pub strike:T,pub notional:T,pub is_call:bool,pub cap:Option<T>}
#[derive(Copy,Clone,Debug)]
pub struct WeatherDerivativeValuation<T: Real>{pub npv:T,pub expected_degree_days:T,pub cumulative_variance:T}
pub fn price_weather_hdd_cdd_option<T: Real>(s:&WeatherDerivativeSpec<T>)->WeatherDerivativeValuation<T>{let w=if s.kind==DegreeDayKind::Hdd{1.0}else{-1.0};let mut mean=T::zero();let mut var=T::zero();let mut prev=0.0;for &t in &s.times{let tt=T::from_f64(t);let seasonal=s.mean_level+s.seasonal_amplitude*((T::from_f64(2.0*std::f64::consts::PI/365.0)*tt+s.seasonal_phase).sin());let mu=seasonal-s.market_price_risk*s.daily_vol/s.mean_reversion.max_val(T::from_f64(1e-8))*(T::one()-(-s.mean_reversion*tt).exp())+(s.initial_temperature-s.mean_level)*(-s.mean_reversion*tt).exp();let vv=s.daily_vol*s.daily_vol/(T::from_f64(2.0)*s.mean_reversion.max_val(T::from_f64(1e-8)))*(T::one()-(-T::from_f64(2.0)*s.mean_reversion*tt).exp());let sd=vv.max_val(T::zero()).sqrt();let z=T::from_f64(w)*(s.temperature_base-mu)/sd.max_val(T::from_f64(1e-8));mean += T::from_f64(w)*(s.temperature_base-mu)*z.norm_cdf()+sd*z.norm_pdf(); var += vv; prev=t;} let st=var.max_val(T::zero()).sqrt(); let d=(mean-s.strike)/st.max_val(T::from_f64(1e-8)); let mut npv=s.notional*(if s.is_call{(mean-s.strike)*d.norm_cdf()+st*d.norm_pdf()}else{(s.strike-mean)*(-d).norm_cdf()+st*d.norm_pdf()}); if let Some(cap)=s.cap{npv=npv.min_val(cap)} let _=prev; WeatherDerivativeValuation{npv,expected_degree_days:mean,cumulative_variance:var}}
#[derive(Clone,Debug)]
pub struct CleanSparkDarkSpreadSpec<T: Real=f64>{pub power:T,pub fuel:T,pub carbon:T,pub heat_rate:T,pub emission_factor:T,pub strike:T,pub vol_power:T,pub vol_fuel:T,pub vol_carbon:T,pub rho_power_fuel:T,pub rho_power_carbon:T,pub rho_fuel_carbon:T,pub rate:T,pub expiry:T,pub option_type:OptionType}
pub fn price_clean_spark_dark_spread<T: Real>(s:&CleanSparkDarkSpreadSpec<T>)->T{if s.emission_factor.val().abs()<1e-14{return spark_spread_option(s.power,s.fuel,s.heat_rate,s.strike,s.vol_power,s.vol_fuel,s.rho_power_fuel,s.rate,s.expiry,s.option_type)} let adj=s.fuel*s.heat_rate+s.carbon*s.emission_factor; kirk_spread_option(s.power,adj,s.strike,s.vol_power,(s.vol_fuel*s.vol_fuel+s.vol_carbon*s.vol_carbon+T::from_f64(2.0)*s.rho_fuel_carbon*s.vol_fuel*s.vol_carbon).sqrt(),s.rho_power_fuel,s.rate,s.expiry,s.option_type)}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    #[test]
    fn clean_spread_zero_carbon_collapse() { let s=CleanSparkDarkSpreadSpec{power:80.0,fuel:3.0,carbon:20.0,heat_rate:7.0,emission_factor:0.0,strike:10.0,vol_power:0.2,vol_fuel:0.3,vol_carbon:0.4,rho_power_fuel:0.2,rho_power_carbon:0.1,rho_fuel_carbon:0.2,rate:0.03,expiry:1.0,option_type:OptionType::Call}; let a=price_clean_spark_dark_spread(&s); let b=spark_spread_option(80.0,3.0,7.0,10.0,0.2,0.3,0.2,0.03,1.0,OptionType::Call); assert!((a-b).abs()<1e-12); }
}
