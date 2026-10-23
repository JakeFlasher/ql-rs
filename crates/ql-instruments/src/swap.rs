use crate::legs::ResolvedCouponLeg;
use ql_ad::Real;
use ql_curves::YieldTermStructure;
use ql_time::{
    BitsetCalendar, BusinessDayConvention, Date, DateGenerationRule, DayCounter, Period,
    ScheduleBuilder,
};
use std::sync::Arc;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SwapType {
    Payer,
    Receiver,
}
#[derive(Debug, Clone)]
pub struct VanillaSwapSpec {
    pub effective_date: Date,
    pub termination_date: Date,
    pub notional: f64,
    pub fixed_rate_or_coupon: f64,
    pub fixed_tenor: Period,
    pub float_tenor: Period,
    pub calendar: Arc<BitsetCalendar>,
    pub bdc: BusinessDayConvention,
    pub fixed_dc: DayCounter,
    pub float_dc: DayCounter,
    pub swap_type: SwapType,
}
#[derive(Debug, Clone)]
pub struct ResolvedVanillaSwap {
    pub fixed_leg: ResolvedCouponLeg,
    pub floating_leg: ResolvedCouponLeg,
    pub swap_type: SwapType,
}
fn leg_from_schedule(
    s: &ql_time::Schedule,
    ref_date: Date,
    notional: f64,
    rate: f64,
    dc: &DayCounter,
) -> ResolvedCouponLeg {
    let mut l = ResolvedCouponLeg::default();
    for i in 1..s.dates.len() {
        let a = s.dates[i - 1];
        let b = s.dates[i];
        let tau = dc.year_fraction(a, b, None, None);
        l.payment_dates.push(b);
        l.accrual_start_dates.push(a);
        l.accrual_end_dates.push(b);
        l.fixing_dates.push(a);
        l.year_fractions.push(tau);
        l.notionals.push(notional);
        l.rates_or_spreads.push(rate);
        l.gearings.push(1.0);
        l.payment_times
            .push(dc.year_fraction(ref_date, b, None, None));
        l.accrual_start_times
            .push(dc.year_fraction(ref_date, a, None, None));
        l.accrual_end_times
            .push(dc.year_fraction(ref_date, b, None, None));
    }
    l.day_counter = Some(dc.clone());
    l
}
impl VanillaSwapSpec {
    pub fn resolve(&self, reference_date: Date) -> ResolvedVanillaSwap {
        let fs = ScheduleBuilder::new(
            self.effective_date,
            self.termination_date,
            self.fixed_tenor,
            self.calendar.clone(),
        )
        .convention(self.bdc)
        .rule(DateGenerationRule::Forward)
        .build();
        let fl = ScheduleBuilder::new(
            self.effective_date,
            self.termination_date,
            self.float_tenor,
            self.calendar.clone(),
        )
        .convention(self.bdc)
        .rule(DateGenerationRule::Forward)
        .build();
        ResolvedVanillaSwap {
            fixed_leg: leg_from_schedule(
                &fs,
                reference_date,
                self.notional,
                self.fixed_rate_or_coupon,
                &self.fixed_dc,
            ),
            floating_leg: leg_from_schedule(
                &fl,
                reference_date,
                self.notional,
                0.0,
                &self.float_dc,
            ),
            swap_type: self.swap_type,
        }
    }
}
impl ResolvedVanillaSwap {
    pub fn npv<T: Real>(
        &self,
        discount: &impl YieldTermStructure<T>,
        forward: &impl YieldTermStructure<T>,
    ) -> T {
        let x = self.fixed_leg.npv_fixed_leg(discount);
        let y = self.floating_leg.npv_floating_leg(discount, forward);
        if matches!(self.swap_type, SwapType::Payer) {
            y - x
        } else {
            x - y
        }
    }
    pub fn annuity<T: Real>(&self, discount: &impl YieldTermStructure<T>) -> T {
        let mut x = T::zero();
        for i in 0..self.fixed_leg.len() {
            x += T::from_f64(self.fixed_leg.notionals[i] * self.fixed_leg.year_fractions[i])
                * discount.discount_factor(self.fixed_leg.payment_times[i]);
        }
        x
    }
    pub fn fair_rate<T: Real>(
        &self,
        discount: &impl YieldTermStructure<T>,
        forward: &impl YieldTermStructure<T>,
    ) -> T {
        self.floating_leg.npv_floating_leg(discount, forward) / self.annuity(discount)
    }
    pub fn fair_spread<T: Real>(
        &self,
        discount: &impl YieldTermStructure<T>,
        forward: &impl YieldTermStructure<T>,
    ) -> T {
        self.fair_rate(discount, forward)
    }
}
#[derive(Debug, Clone)]
pub struct OisSwapSpec {
    pub effective_date: Date,
    pub termination_date: Date,
    pub notional: f64,
    pub fixed_rate_or_coupon: f64,
    pub tenor: Period,
    pub calendar: Arc<BitsetCalendar>,
    pub bdc: BusinessDayConvention,
    pub fixed_dc: DayCounter,
    pub float_dc: DayCounter,
    pub swap_type: SwapType,
}
#[derive(Debug, Clone)]
pub struct ResolvedOisSwap {
    pub inner: ResolvedVanillaSwap,
}
impl OisSwapSpec {
    pub fn resolve(&self, reference_date: Date) -> ResolvedOisSwap {
        let v = VanillaSwapSpec {
            effective_date: self.effective_date,
            termination_date: self.termination_date,
            notional: self.notional,
            fixed_rate_or_coupon: self.fixed_rate_or_coupon,
            fixed_tenor: self.tenor,
            float_tenor: self.tenor,
            calendar: self.calendar.clone(),
            bdc: self.bdc,
            fixed_dc: self.fixed_dc.clone(),
            float_dc: self.float_dc.clone(),
            swap_type: self.swap_type,
        };
        ResolvedOisSwap {
            inner: v.resolve(reference_date),
        }
    }
}
impl ResolvedOisSwap {
    pub fn npv<T: Real>(&self, curve: &impl YieldTermStructure<T>) -> T {
        self.inner.npv(curve, curve)
    }
    pub fn fair_rate<T: Real>(&self, curve: &impl YieldTermStructure<T>) -> T {
        self.inner.fair_rate(curve, curve)
    }
}
#[derive(Debug, Clone)]
pub struct BasisSwapSpec {
    pub effective_date: Date,
    pub termination_date: Date,
    pub notional: f64,
    pub tenor1: Period,
    pub tenor2: Period,
    pub spread: f64,
    pub calendar: Arc<BitsetCalendar>,
    pub bdc: BusinessDayConvention,
    pub day_counter: DayCounter,
    pub swap_type: SwapType,
}
#[derive(Debug, Clone)]
pub struct ResolvedBasisSwap {
    pub leg1: ResolvedCouponLeg,
    pub leg2: ResolvedCouponLeg,
    pub swap_type: SwapType,
}
impl BasisSwapSpec {
    pub fn resolve(&self, reference_date: Date) -> ResolvedBasisSwap {
        let s1 = ScheduleBuilder::new(
            self.effective_date,
            self.termination_date,
            self.tenor1,
            self.calendar.clone(),
        )
        .convention(self.bdc)
        .rule(DateGenerationRule::Forward)
        .build();
        let s2 = ScheduleBuilder::new(
            self.effective_date,
            self.termination_date,
            self.tenor2,
            self.calendar.clone(),
        )
        .convention(self.bdc)
        .rule(DateGenerationRule::Forward)
        .build();
        ResolvedBasisSwap {
            leg1: leg_from_schedule(
                &s1,
                reference_date,
                self.notional,
                self.spread,
                &self.day_counter,
            ),
            leg2: leg_from_schedule(&s2, reference_date, self.notional, 0.0, &self.day_counter),
            swap_type: self.swap_type,
        }
    }
}
impl ResolvedBasisSwap {
    pub fn npv<T: Real>(
        &self,
        discount: &impl YieldTermStructure<T>,
        fwd1: &impl YieldTermStructure<T>,
        fwd2: &impl YieldTermStructure<T>,
    ) -> T {
        let a = self.leg1.npv_floating_leg(discount, fwd1);
        let b = self.leg2.npv_floating_leg(discount, fwd2);
        if matches!(self.swap_type, SwapType::Payer) {
            a - b
        } else {
            b - a
        }
    }
}
