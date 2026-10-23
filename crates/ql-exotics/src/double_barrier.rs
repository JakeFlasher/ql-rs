use ql_ad::Real;
use ql_instruments::OptionType;
use ql_vol::black_scholes_price;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum DoubleBarrierType {
    KnockOut,
    KnockIn,
}

#[derive(Copy, Clone, Debug)]
pub struct DoubleBarrierOptionSpec<T: Real> {
    pub option_type: OptionType,
    pub barrier_type: DoubleBarrierType,
    pub spot: T,
    pub strike: T,
    pub barrier_lo: T,
    pub barrier_hi: T,
    pub rate: T,
    pub div_yield: T,
    pub vol: T,
    pub expiry: T,
    pub series_terms: i32,
}

pub fn price_double_barrier_option<T: Real>(spec: &DoubleBarrierOptionSpec<T>) -> T {
    let vanilla = black_scholes_price(
        spec.spot,
        spec.strike,
        spec.vol,
        spec.expiry,
        spec.rate,
        spec.div_yield,
        spec.option_type.is_call(),
    );
    if spec.spot.val() <= spec.barrier_lo.val() || spec.spot.val() >= spec.barrier_hi.val() {
        return match spec.barrier_type {
            DoubleBarrierType::KnockOut => T::zero(),
            DoubleBarrierType::KnockIn => vanilla,
        };
    }
    if spec.expiry.val() <= 0.0 || spec.vol.val() <= 0.0 {
        return match spec.barrier_type {
            DoubleBarrierType::KnockOut => T::zero(),
            DoubleBarrierType::KnockIn => vanilla,
        };
    }
    let std = spec.vol * spec.expiry.sqrt();
    let sig2 = spec.vol * spec.vol;
    let b = spec.rate - spec.div_yield;
    let mu1 = b * 2.0 / sig2 + 1.0;
    let bsig = (b + sig2 * 0.5) * spec.expiry / std;
    let ln_ratio = (spec.barrier_hi / spec.barrier_lo).ln();
    let ln_l = spec.barrier_lo.ln();
    let ln_u = spec.barrier_hi.ln();
    let ln_s = spec.spot.ln();
    let ns = spec.series_terms.max(1);
    let mut acc1 = T::zero();
    let mut acc2 = T::zero();
    let add = |acc: &mut T, w: T, p: T, sign: f64| {
        if p.val().abs() > 0.0 {
            *acc += T::from_f64(sign) * w * p;
        }
    };
    if spec.option_type.is_call() {
        if spec.strike.val() >= spec.barrier_hi.val() {
            acc1 = T::zero();
            acc2 = T::zero();
        } else {
            let bound = spec.strike.max_val(spec.barrier_lo);
            for n in -ns..=ns {
                let nf = n as f64;
                let shift = ln_ratio * (2.0 * nf);
                let d1 = ((spec.spot / bound).ln() + shift) / std + bsig;
                let d2 = ((spec.spot / spec.barrier_hi).ln() + shift) / std + bsig;
                let d3 = ((spec.barrier_lo * spec.barrier_lo / (bound * spec.spot)).ln() - shift)
                    / std
                    + bsig;
                let d4 = ((spec.barrier_lo * spec.barrier_lo / (spec.barrier_hi * spec.spot)).ln()
                    - shift)
                    / std
                    + bsig;
                let w1a = (mu1 * nf * ln_ratio).exp();
                let w1b = (mu1 * (ln_l * (nf + 1.0) - ln_u * nf - ln_s)).exp();
                let mu2 = mu1 - 2.0;
                let w2a = (mu2 * nf * ln_ratio).exp();
                let w2b = (mu2 * (ln_l * (nf + 1.0) - ln_u * nf - ln_s)).exp();
                add(&mut acc1, w1a, d1.norm_cdf() - d2.norm_cdf(), 1.0);
                add(&mut acc1, w1b, d3.norm_cdf() - d4.norm_cdf(), -1.0);
                add(
                    &mut acc2,
                    w2a,
                    (d1 - std).norm_cdf() - (d2 - std).norm_cdf(),
                    1.0,
                );
                add(
                    &mut acc2,
                    w2b,
                    (d3 - std).norm_cdf() - (d4 - std).norm_cdf(),
                    -1.0,
                );
            }
        }
    } else if spec.strike.val() <= spec.barrier_lo.val() {
        acc1 = T::zero();
        acc2 = T::zero();
    } else {
        let bound = spec.strike.min_val(spec.barrier_hi);
        for n in -ns..=ns {
            let nf = n as f64;
            let shift = ln_ratio * (2.0 * nf);
            let y1 = ((spec.spot / spec.barrier_lo).ln() + shift) / std + bsig;
            let y2 = ((spec.spot / bound).ln() + shift) / std + bsig;
            let y3 = ((spec.barrier_lo / spec.spot).ln() - shift) / std + bsig;
            let y4 = ((spec.barrier_lo * spec.barrier_lo / (bound * spec.spot)).ln() - shift) / std
                + bsig;
            let w1a = (mu1 * nf * ln_ratio).exp();
            let w1b = (mu1 * (ln_l * (nf + 1.0) - ln_u * nf - ln_s)).exp();
            let mu2 = mu1 - 2.0;
            let w2a = (mu2 * nf * ln_ratio).exp();
            let w2b = (mu2 * (ln_l * (nf + 1.0) - ln_u * nf - ln_s)).exp();
            add(
                &mut acc1,
                w2a,
                (y1 - std).norm_cdf() - (y2 - std).norm_cdf(),
                1.0,
            );
            add(
                &mut acc1,
                w2b,
                (y3 - std).norm_cdf() - (y4 - std).norm_cdf(),
                -1.0,
            );
            add(&mut acc2, w1a, y1.norm_cdf() - y2.norm_cdf(), 1.0);
            add(&mut acc2, w1b, y3.norm_cdf() - y4.norm_cdf(), -1.0);
        }
    }
    let ko = if spec.option_type.is_call() {
        (spec.spot * (-spec.div_yield * spec.expiry).exp() * acc1
            - spec.strike * (-spec.rate * spec.expiry).exp() * acc2)
            .max_val(T::zero())
    } else {
        (spec.strike * (-spec.rate * spec.expiry).exp() * acc1
            - spec.spot * (-spec.div_yield * spec.expiry).exp() * acc2)
            .max_val(T::zero())
    };
    match spec.barrier_type {
        DoubleBarrierType::KnockOut => ko,
        DoubleBarrierType::KnockIn => vanilla - ko,
    }
}
