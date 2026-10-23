use ql_ad::{Dual, Real};
use ql_instruments::OptionType;
use ql_vol::black_scholes_price;

pub const BGK_BETA: f64 = 0.5825971579390107;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BarrierType {
    DownIn,
    UpIn,
    DownOut,
    UpOut,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MonitoringFrequency {
    Continuous,
    Discrete { observations: usize },
}

#[derive(Clone, Copy, Debug)]
pub struct BarrierOptionSpec<T: Real> {
    pub option_type: OptionType,
    pub barrier_type: BarrierType,
    pub spot: T,
    pub strike: T,
    pub barrier: T,
    pub rebate: T,
    pub rate: T,
    pub div_yield: T,
    pub vol: T,
    pub expiry: T,
    pub monitoring: MonitoringFrequency,
}

#[inline]
pub fn bgk_adjusted_barrier<T: Real>(
    spot: T,
    barrier: T,
    vol: T,
    expiry: T,
    observations: usize,
) -> T {
    if observations == 0 {
        return barrier;
    }
    let dt = expiry / observations as f64;
    let shift = vol * dt.sqrt() * BGK_BETA;
    if barrier.val() > spot.val() {
        barrier * shift.exp()
    } else if barrier.val() < spot.val() {
        barrier * (-shift).exp()
    } else {
        barrier
    }
}

fn price_continuous<T: Real>(spec: &BarrierOptionSpec<T>, h: T) -> T {
    let s = spec.spot;
    let k = spec.strike;
    let r = spec.rate;
    let q = spec.div_yield;
    let v = spec.vol;
    let t = spec.expiry;
    if t.val() <= 0.0 || v.val() <= 0.0 {
        return black_scholes_price(s, k, v, t, r, q, spec.option_type.is_call());
    }
    let b = r - q;
    let vv = v * v;
    let st = v * t.sqrt();
    let mu = (b - T::from_f64(0.5) * vv) / vv;
    let lambda = (mu * mu + T::from_f64(2.0) * r / vv).sqrt();
    let x1 = (s / k).ln() / st + (T::one() + mu) * st;
    let x2 = (s / h).ln() / st + (T::one() + mu) * st;
    let y1 = (h * h / (s * k)).ln() / st + (T::one() + mu) * st;
    let y2 = (h / s).ln() / st + (T::one() + mu) * st;
    let z = (h / s).ln() / st + lambda * st;
    let dr = (-r * t).exp();
    let dq = ((b - r) * t).exp();
    let ratio = h / s;
    let phi = if spec.option_type.is_call() {
        1.0
    } else {
        -1.0
    };
    let eta = if matches!(
        spec.barrier_type,
        BarrierType::DownIn | BarrierType::DownOut
    ) {
        1.0
    } else {
        -1.0
    };
    let a = T::from_f64(phi)
        * (s * dq * (T::from_f64(phi) * x1).norm_cdf()
            - k * dr * (T::from_f64(phi) * (x1 - st)).norm_cdf());
    let bb = T::from_f64(phi)
        * (s * dq * (T::from_f64(phi) * x2).norm_cdf()
            - k * dr * (T::from_f64(phi) * (x2 - st)).norm_cdf());
    let c = T::from_f64(phi)
        * (s * dq * ratio.pow((mu + 1.0) * 2.0) * (T::from_f64(eta) * y1).norm_cdf()
            - k * dr * ratio.pow(mu * 2.0) * (T::from_f64(eta) * (y1 - st)).norm_cdf());
    let d = T::from_f64(phi)
        * (s * dq * ratio.pow((mu + 1.0) * 2.0) * (T::from_f64(eta) * y2).norm_cdf()
            - k * dr * ratio.pow(mu * 2.0) * (T::from_f64(eta) * (y2 - st)).norm_cdf());
    let e = spec.rebate
        * dr
        * ((T::from_f64(eta) * (x2 - st)).norm_cdf()
            - ratio.pow(mu * 2.0) * (T::from_f64(eta) * (y2 - st)).norm_cdf());
    let f = spec.rebate
        * (ratio.pow(mu + lambda) * (T::from_f64(eta) * z).norm_cdf()
            + ratio.pow(mu - lambda)
                * (T::from_f64(eta) * (z - T::from_f64(2.0) * lambda * st)).norm_cdf());
    let down = matches!(
        spec.barrier_type,
        BarrierType::DownIn | BarrierType::DownOut
    );
    let in_barrier = matches!(spec.barrier_type, BarrierType::DownIn | BarrierType::UpIn);
    let ge = k.val() >= h.val();
    match (down, in_barrier, spec.option_type.is_call(), ge) {
        (true, true, true, true) => c + e,
        (true, true, true, false) => a - bb + d + e,
        (true, false, true, true) => a - c + f,
        (true, false, true, false) => bb - d + f,
        (false, true, true, true) => a + e,
        (false, true, true, false) => bb - c + d + e,
        (false, false, true, true) => f,
        (false, false, true, false) => a - bb + c - d + f,
        (true, true, false, true) => bb - c + d + e,
        (true, true, false, false) => a + e,
        (true, false, false, true) => a - bb + c - d + f,
        (true, false, false, false) => f,
        (false, true, false, true) => a - bb + d + e,
        (false, true, false, false) => c + e,
        (false, false, false, true) => bb - d + f,
        (false, false, false, false) => a - c + f,
    }
}

pub fn price_barrier_option<T: Real>(spec: &BarrierOptionSpec<T>) -> T {
    let adjusted = match spec.monitoring {
        MonitoringFrequency::Continuous => spec.barrier,
        MonitoringFrequency::Discrete { observations } => {
            bgk_adjusted_barrier(spec.spot, spec.barrier, spec.vol, spec.expiry, observations)
        }
    };
    let down = matches!(
        spec.barrier_type,
        BarrierType::DownIn | BarrierType::DownOut
    );
    let breached = if down {
        spec.spot.val() <= adjusted.val()
    } else {
        spec.spot.val() >= adjusted.val()
    };
    let knocked_in = matches!(spec.barrier_type, BarrierType::DownIn | BarrierType::UpIn);
    if breached {
        return if knocked_in {
            black_scholes_price(
                spec.spot,
                spec.strike,
                spec.vol,
                spec.expiry,
                spec.rate,
                spec.div_yield,
                spec.option_type.is_call(),
            )
        } else {
            spec.rebate
        };
    }
    price_continuous(spec, adjusted)
}

pub fn barrier_option_greeks_dual5(spec: &BarrierOptionSpec<f64>) -> (f64, [f64; 5]) {
    let d = BarrierOptionSpec {
        option_type: spec.option_type,
        barrier_type: spec.barrier_type,
        spot: Dual::<5>::variable(spec.spot, 0),
        strike: Dual::from_value(spec.strike),
        barrier: Dual::from_value(spec.barrier),
        rebate: Dual::from_value(spec.rebate),
        rate: Dual::variable(spec.rate, 2),
        div_yield: Dual::variable(spec.div_yield, 3),
        vol: Dual::variable(spec.vol, 1),
        expiry: Dual::variable(spec.expiry, 4),
        monitoring: spec.monitoring,
    };
    let p = price_barrier_option(&d);
    (p.val, [p.d[0], p.d[1], p.d[2], p.d[3], -p.d[4]])
}
