use ql_ad::{Real, SimdReal};
use ql_curves::YieldTermStructure;
use ql_time::{Date, DayCounter};
#[derive(Debug, Clone, Default)]
pub struct ResolvedCouponLeg {
    pub payment_dates: Vec<Date>,
    pub accrual_start_dates: Vec<Date>,
    pub accrual_end_dates: Vec<Date>,
    pub fixing_dates: Vec<Date>,
    pub year_fractions: Vec<f64>,
    pub notionals: Vec<f64>,
    pub rates_or_spreads: Vec<f64>,
    pub gearings: Vec<f64>,
    pub payment_times: Vec<f64>,
    pub accrual_start_times: Vec<f64>,
    pub accrual_end_times: Vec<f64>,
    pub day_counter: Option<DayCounter>,
}

/// Compact fixed-rate SoA view used by high-throughput flat-curve valuation.
/// Coupon amounts are materialized once at resolution, so a scenario only
/// loads contiguous cashflows and payment times.
#[derive(Debug, Clone)]
pub struct FixedRateLeg {
    pub cashflows: Vec<f64>,
    pub payment_times: Vec<f64>,
    pub regular_dt: Option<f64>,
}

impl FixedRateLeg {
    pub fn from_resolved(leg: &ResolvedCouponLeg) -> Self {
        assert_eq!(leg.notionals.len(), leg.year_fractions.len());
        assert_eq!(leg.notionals.len(), leg.rates_or_spreads.len());
        assert_eq!(leg.notionals.len(), leg.payment_times.len());
        let cashflows = leg
            .notionals
            .iter()
            .zip(&leg.year_fractions)
            .zip(&leg.rates_or_spreads)
            .map(|((n, yf), r)| n * yf * r)
            .collect::<Vec<_>>();
        let regular_dt = if leg.payment_times.len() > 1 {
            let dt = leg.payment_times[1] - leg.payment_times[0];
            (dt > 0.0
                && leg
                    .payment_times
                    .windows(2)
                    .all(|w| (w[1] - w[0] - dt).abs() < 1e-14))
            .then_some(dt)
        } else {
            None
        };
        Self {
            cashflows,
            payment_times: leg.payment_times.clone(),
            regular_dt,
        }
    }

    #[inline(always)]
    pub fn npv_simd_flat(&self, rate: f64) -> f64 {
        use core::simd::Simd;
        use core::simd::prelude::SimdFloat;
        let neg_rate = Simd::<f64, 8>::splat(-rate);
        let mut total = Simd::<f64, 8>::splat(0.0);
        let mut i = 0;
        while i + 8 <= self.cashflows.len() {
            let cash = Simd::<f64, 8>::from_slice(&self.cashflows[i..i + 8]);
            let times = Simd::<f64, 8>::from_slice(&self.payment_times[i..i + 8]);
            total += cash * SimdReal::exp(times * neg_rate);
            i += 8;
        }
        let mut out = total.reduce_sum();
        while i < self.cashflows.len() {
            out += self.cashflows[i] * (-rate * self.payment_times[i]).exp();
            i += 1;
        }
        out
    }

    /// Value independent flat-rate scenarios in eight SIMD lanes.
    pub fn npv_simd_flat_batch(&self, rates: &[f64]) -> Vec<f64> {
        use core::simd::Simd;
        let mut out = vec![0.0; rates.len()];
        for (batch, chunk) in rates.chunks(8).enumerate() {
            let mut rate_arr = [0.0; 8];
            for (lane, rate) in chunk.iter().enumerate() {
                rate_arr[lane] = *rate;
            }
            let neg_rate = -Simd::<f64, 8>::from_array(rate_arr);
            let mut sum = Simd::<f64, 8>::splat(0.0);
            if let Some(dt) = self.regular_dt {
                let mut discount = SimdReal::exp(neg_rate * Simd::splat(self.payment_times[0]));
                let q = SimdReal::exp(neg_rate * Simd::splat(dt));
                for &cash in &self.cashflows {
                    sum += discount * Simd::splat(cash);
                    discount *= q;
                }
            } else {
                for (&cash, &time) in self.cashflows.iter().zip(&self.payment_times) {
                    sum += SimdReal::exp(neg_rate * Simd::splat(time)) * Simd::splat(cash);
                }
            }
            let values = sum.to_array();
            let start = batch * 8;
            out[start..start + chunk.len()].copy_from_slice(&values[..chunk.len()]);
        }
        out
    }
}

impl ResolvedCouponLeg {
    pub fn is_empty(&self) -> bool {
        self.payment_dates.is_empty()
    }
    pub fn len(&self) -> usize {
        self.payment_dates.len()
    }
    pub fn npv_fixed_leg<T: Real>(&self, discount: &impl YieldTermStructure<T>) -> T {
        let mut x = T::zero();
        for i in 0..self.len() {
            let t = self.payment_times.get(i).copied().unwrap_or(0.0);
            x += T::from_f64(self.notionals[i] * self.year_fractions[i] * self.rates_or_spreads[i])
                * discount.discount_factor(t);
        }
        x
    }
    pub fn npv_floating_leg<T: Real>(
        &self,
        discount: &impl YieldTermStructure<T>,
        forward: &impl YieldTermStructure<T>,
    ) -> T {
        let mut x = T::zero();
        for i in 0..self.len() {
            let a = self.accrual_start_times.get(i).copied().unwrap_or(0.0);
            let b = self
                .accrual_end_times
                .get(i)
                .copied()
                .unwrap_or(a + self.year_fractions[i]);
            let f = (forward.discount_factor(a) / forward.discount_factor(b) - T::one())
                / self.year_fractions[i];
            let rate = f * self.gearings.get(i).copied().unwrap_or(1.0) + self.rates_or_spreads[i];
            x += T::from_f64(self.notionals[i] * self.year_fractions[i])
                * rate
                * discount.discount_factor(self.payment_times.get(i).copied().unwrap_or(b));
        }
        x
    }
    pub fn npv_fixed_leg_simd(&self, discount: &impl YieldTermStructure<f64>) -> f64 {
        use core::simd::Simd;
        use core::simd::prelude::SimdFloat;
        let mut total = Simd::<f64, 8>::splat(0.0);
        let mut i = 0;
        while i + 8 <= self.len() {
            let notionals = Simd::<f64, 8>::from_slice(&self.notionals[i..i + 8]);
            let year_fractions = Simd::<f64, 8>::from_slice(&self.year_fractions[i..i + 8]);
            let rates = Simd::<f64, 8>::from_slice(&self.rates_or_spreads[i..i + 8]);
            let discounts = Simd::<f64, 8>::from_array(std::array::from_fn(|lane| {
                discount.discount_factor(self.payment_times[i + lane])
            }));
            total += notionals * year_fractions * rates * discounts;
            i += 8;
        }
        let mut out = total.reduce_sum();
        while i < self.len() {
            out += self.notionals[i]
                * self.year_fractions[i]
                * self.rates_or_spreads[i]
                * discount.discount_factor(self.payment_times[i]);
            i += 1;
        }
        out
    }

    /// SIMD valuation against a flat continuously compounded discount rate.
    pub fn npv_simd_flat(&self, rate: f64) -> f64 {
        use core::simd::Simd;
        use core::simd::prelude::SimdFloat;
        let mut total = Simd::<f64, 8>::splat(0.0);
        let neg_rate = Simd::<f64, 8>::splat(-rate);
        let mut i = 0;
        while i + 8 <= self.len() {
            let notionals = Simd::<f64, 8>::from_slice(&self.notionals[i..i + 8]);
            let year_fractions = Simd::<f64, 8>::from_slice(&self.year_fractions[i..i + 8]);
            let rates = Simd::<f64, 8>::from_slice(&self.rates_or_spreads[i..i + 8]);
            let times = Simd::<f64, 8>::from_slice(&self.payment_times[i..i + 8]);
            let cash = notionals * year_fractions * rates;
            total += cash * SimdReal::exp(times * neg_rate);
            i += 8;
        }
        let mut out = total.reduce_sum();
        while i < self.len() {
            out += self.notionals[i]
                * self.year_fractions[i]
                * self.rates_or_spreads[i]
                * (-rate * self.payment_times[i]).exp();
            i += 1;
        }
        out
    }
    pub fn npv_swap_fast_f64(&self, discount: &impl YieldTermStructure<f64>) -> f64 {
        self.npv_fixed_leg(discount)
    }
}
impl ResolvedCouponLeg {
    pub fn single(
        payment_date: Date,
        payment_time: f64,
        notional: f64,
        rate: f64,
        year_fraction: f64,
    ) -> Self {
        Self {
            payment_dates: vec![payment_date],
            accrual_start_dates: vec![payment_date],
            accrual_end_dates: vec![payment_date],
            fixing_dates: vec![payment_date],
            year_fractions: vec![year_fraction],
            notionals: vec![notional],
            rates_or_spreads: vec![rate],
            gearings: vec![1.0],
            payment_times: vec![payment_time],
            accrual_start_times: vec![payment_time - year_fraction],
            accrual_end_times: vec![payment_time],
            day_counter: None,
        }
    }
}
pub fn npv_swap_fast_f64(leg: &ResolvedCouponLeg, discount: &impl YieldTermStructure<f64>) -> f64 {
    leg.npv_fixed_leg(discount)
}
