use crate::{BitsetCalendar, Date, Frequency, days_in_month};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub enum DayCounter {
    Act360,
    #[default]
    Act365Fixed,
    Act365Leap,
    ActActIsda,
    ActActIcma {
        frequency: Frequency,
    },
    Thirty360Us,
    Thirty360Euro,
    Thirty360Isda,
    #[serde(skip)]
    Business252(Arc<BitsetCalendar>),
}
impl DayCounter {
    pub fn day_count(&self, d1: Date, d2: Date) -> i32 {
        match self {
            Self::Business252(c) => c.business_days_between(d1, d2),
            _ => d2.days().saturating_sub(d1.days()),
        }
    }
    pub fn year_fraction(
        &self,
        d1: Date,
        d2: Date,
        ref_start: Option<Date>,
        ref_end: Option<Date>,
    ) -> f64 {
        if d1 == d2 {
            return 0.0;
        }
        let sign = if d2 > d1 { 1.0 } else { -1.0 };
        let (a, b) = if sign > 0.0 { (d1, d2) } else { (d2, d1) };
        let value = match self {
            Self::Act360 => (i64::from(b.days()) - i64::from(a.days())) as f64 / 360.0,
            Self::Act365Fixed => (i64::from(b.days()) - i64::from(a.days())) as f64 / 365.0,
            Self::Act365Leap => act365_l(a, b),
            Self::ActActIsda => act365_isda(a, b),
            Self::ActActIcma { frequency } => {
                let rs = ref_start.unwrap_or(a);
                let re = ref_end
                    .unwrap_or_else(|| rs.add_months(12 / frequency.periods_per_year(), false));
                (i64::from(b.days()) - i64::from(a.days())) as f64
                    / ((i64::from(re.days()) - i64::from(rs.days())) as f64
                        * frequency.periods_per_year() as f64)
            }
            Self::Thirty360Us => thirty_us(a, b) as f64 / 360.0,
            Self::Thirty360Euro => thirty_euro(a, b) as f64 / 360.0,
            Self::Thirty360Isda => thirty_isda(a, b) as f64 / 360.0,
            Self::Business252(c) => c.business_days_between(a, b) as f64 / 252.0,
        };
        sign * value
    }
}
fn act365_l(a: Date, b: Date) -> f64 {
    let days = (i64::from(b.days()) - i64::from(a.days())) as f64;
    let mut y = a.year();
    while y <= b.year() {
        if Date::is_leap_year(y)
            && let Ok(feb29) = Date::from_ymd(y, 2, 29)
            && feb29 >= a
            && feb29 < b
        {
            return days / 366.0;
        }
        y += 1;
    }
    days / 365.0
}
fn act365_isda(a: Date, b: Date) -> f64 {
    let mut total = 0.0;
    let mut y = a.year();
    while y <= b.year() {
        let ys = Date::from_ymd(y, 1, 1).expect("valid year");
        let ye = Date::from_ymd(y + 1, 1, 1).expect("valid year");
        let lo = a.max(ys);
        let hi = b.min(ye);
        if hi > lo {
            total += (i64::from(hi.days()) - i64::from(lo.days())) as f64
                / if Date::is_leap_year(y) { 366.0 } else { 365.0 };
        }
        y += 1;
    }
    total
}
fn thirty_us(a: Date, b: Date) -> i32 {
    let (_, am, mut ad) = a.ymd();
    let (_, bm, mut bd) = b.ymd();
    let a_last_feb = am == 2 && ad == days_in_month(a.year(), am);
    let b_last_feb = bm == 2 && bd == days_in_month(b.year(), bm);
    if ad == 31 || a_last_feb {
        ad = 30;
    }
    // Under the NASD 30/360 rule, an end-of-February date is treated as day
    // 30 only when the start is also end-of-February (which has already been
    // normalised above).  This preserves the special February-to-February
    // convention while leaving a March/February end date untouched.
    if (bd == 31 && ad == 30) || (b_last_feb && a_last_feb) {
        bd = 30;
    }
    ((i64::from(b.year()) - i64::from(a.year())) * 360
        + (i64::from(bm) - i64::from(am)) * 30
        + i64::from(bd)
        - i64::from(ad))
    .clamp(i32::MIN as i64, i32::MAX as i64) as i32
}
fn thirty_euro(a: Date, b: Date) -> i32 {
    let (_, am, mut ad) = a.ymd();
    let (_, bm, mut bd) = b.ymd();
    if ad == 31 {
        ad = 30;
    }
    if bd == 31 {
        bd = 30;
    }
    ((i64::from(b.year()) - i64::from(a.year())) * 360
        + (i64::from(bm) - i64::from(am)) * 30
        + i64::from(bd)
        - i64::from(ad))
    .clamp(i32::MIN as i64, i32::MAX as i64) as i32
}
fn thirty_isda(a: Date, b: Date) -> i32 {
    let (_, am, mut ad) = a.ymd();
    let (_, bm, mut bd) = b.ymd();
    if ad == 31 || (am == 2 && ad == days_in_month(a.year(), am)) {
        ad = 30;
    }
    if bd == 31 || (bm == 2 && bd == days_in_month(b.year(), bm)) {
        bd = 30;
    }
    ((i64::from(b.year()) - i64::from(a.year())) * 360
        + (i64::from(bm) - i64::from(am)) * 30
        + i64::from(bd)
        - i64::from(ad))
    .clamp(i32::MIN as i64, i32::MAX as i64) as i32
}
impl std::ops::Sub for Date {
    type Output = i32;
    fn sub(self, rhs: Self) -> i32 {
        self.days().saturating_sub(rhs.days())
    }
}
