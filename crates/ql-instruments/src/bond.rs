use ql_ad::{Dual2, Real, ift_lift_1d};
use ql_math::newton_solve;
use ql_time::{
    BitsetCalendar, BusinessDayConvention, Date, DateGenerationRule, DayCounter, Frequency, Period,
    ScheduleBuilder,
};
use std::sync::Arc;
#[derive(Debug, Clone)]
pub struct FixedRateBondSpec {
    pub effective_date: Date,
    pub maturity_date: Date,
    pub coupon_rate: f64,
    pub frequency: Frequency,
    pub notional: f64,
    pub calendar: Arc<BitsetCalendar>,
    pub bdc: BusinessDayConvention,
    pub day_counter: DayCounter,
}
#[derive(Debug, Clone)]
pub struct ResolvedFixedRateBond {
    pub settlement_date: Date,
    pub maturity_date: Date,
    pub coupon_rate: f64,
    pub frequency: Frequency,
    pub notional: f64,
    pub payment_times: Vec<f64>,
    pub payment_dates: Vec<Date>,
    pub accrual_start_dates: Vec<Date>,
    pub accrual_end_dates: Vec<Date>,
    pub day_counter: DayCounter,
}
impl FixedRateBondSpec {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        effective_date: Date,
        maturity_date: Date,
        coupon_rate: f64,
        frequency: Frequency,
        notional: f64,
        calendar: Arc<BitsetCalendar>,
        bdc: BusinessDayConvention,
        day_counter: DayCounter,
    ) -> Self {
        Self {
            effective_date,
            maturity_date,
            coupon_rate,
            frequency,
            notional,
            calendar,
            bdc,
            day_counter,
        }
    }
    pub fn resolve(&self, reference_date: Date) -> ResolvedFixedRateBond {
        let s = ScheduleBuilder::new(
            self.effective_date,
            self.maturity_date,
            Period::from(self.frequency),
            self.calendar.clone(),
        )
        .convention(self.bdc)
        .rule(DateGenerationRule::Forward)
        .build();
        let mut p = Vec::new();
        let mut a = Vec::new();
        let mut b = Vec::new();
        for i in 1..s.dates.len() {
            a.push(s.dates[i - 1]);
            b.push(s.dates[i]);
            p.push(self.day_counter.year_fraction(reference_date, s.dates[i], None, None));
        }
        ResolvedFixedRateBond {
            settlement_date: reference_date,
            maturity_date: self.maturity_date,
            coupon_rate: self.coupon_rate,
            frequency: self.frequency,
            notional: self.notional,
            payment_times: p,
            payment_dates: s.dates[1..].to_vec(),
            accrual_start_dates: a,
            accrual_end_dates: b,
            day_counter: self.day_counter.clone(),
        }
    }
}
impl ResolvedFixedRateBond {
    fn freq(&self) -> f64 {
        self.frequency.periods_per_year() as f64
    }
    fn flows<T: Real>(&self, y: T) -> T {
        let c = T::from_f64(self.notional * self.coupon_rate / self.freq());
        let mut p = T::zero();
        for (i, t) in self.payment_times.iter().enumerate() {
            let cf = if i + 1 == self.payment_times.len() {
                c + T::from_f64(self.notional)
            } else {
                c
            };
            p += cf / (T::one() + y / self.freq()).powf(*t * self.freq());
        }
        p
    }
    pub fn dirty_price<T: Real>(&self, y: T) -> T {
        self.flows(y)
    }
    pub fn accrued_interest(&self, settlement: Date) -> f64 {
        for i in 0..self.accrual_end_dates.len() {
            if settlement >= self.accrual_start_dates[i] && settlement < self.accrual_end_dates[i] {
                let tau = self.day_counter.year_fraction(
                    self.accrual_start_dates[i],
                    settlement,
                    None,
                    None,
                );
                return self.notional * self.coupon_rate / self.freq() * tau
                    / self.day_counter.year_fraction(
                        self.accrual_start_dates[i],
                        self.accrual_end_dates[i],
                        None,
                        None,
                    );
            }
        }
        0.0
    }
    pub fn clean_price<T: Real>(&self, y: T) -> T {
        self.dirty_price(y) - T::from_f64(self.accrued_interest(self.settlement_date))
    }
    pub fn yield_to_maturity<T: Real>(&self, dirty_price: T) -> T {
        let p = dirty_price.val();
        let guess = self.coupon_rate.max(-0.99);
        let root = match newton_solve(
            |y| {
                let yy = Dual2::<1>::variable(y, 0);
                let flow = self.flows(yy);
                (flow.val - p, flow.d1[0])
            },
            guess,
            1e-13,
            100,
        ) {
            Ok(x) => x,
            Err(_) => {
                let mut x = guess;
                for _ in 0..64 {
                    let yy = Dual2::<1>::variable(x, 0);
                    let flow = self.flows(yy);
                    if flow.d1[0].abs() < 1e-14 || !flow.d1[0].is_finite() {
                        break;
                    }
                    let xn = x - (flow.val - p) / flow.d1[0];
                    if !xn.is_finite() || (xn - x).abs() < 1e-13 {
                        break;
                    }
                    x = xn;
                }
                x
            }
        };
        let root_ad = Dual2::<1>::variable(root, 0);
        let root_flow = self.flows(root_ad);
        ift_lift_1d(
            root,
            |x: T| self.flows(x) - dirty_price,
            T::from_f64(root_flow.d1[0]),
            root_flow.d2[0][0],
        )
    }
    pub fn modified_duration(&self, y: f64) -> f64 {
        let x = Dual2::<1>::variable(y, 0);
        let p = self.flows(x);
        -p.d1[0] / p.val
    }
    pub fn macaulay_duration(&self, y: f64) -> f64 {
        self.modified_duration(y) * (1.0 + y / self.freq())
    }
    pub fn convexity(&self, y: f64) -> f64 {
        let x = Dual2::<1>::variable(y, 0);
        let p = self.flows(x);
        p.d2[0][0] / p.val
    }
}
