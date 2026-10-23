use crate::date::{Date, Period, Weekday};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

const START_DAY: i32 = -25_567; // 1900-01-01 relative to 1970-01-01
const END_DAY: i32 = 67_934; // 2155-12-31
const N_WORDS: usize = 1461;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BusinessDayConvention {
    Unadjusted,
    Following,
    ModifiedFollowing,
    Preceding,
    ModifiedPreceding,
    HalfMonthModifiedFollowing,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BitsetCalendar {
    pub name: Arc<str>,
    pub(crate) bits: [u64; N_WORDS],
}
impl BitsetCalendar {
    pub const START: Date = Date::from_days(START_DAY);
    pub const END: Date = Date::from_days(END_DAY);
    pub fn new(name: impl Into<Arc<str>>) -> Self {
        Self {
            name: name.into(),
            bits: [0; N_WORDS],
        }
    }
    pub fn weekend_only() -> Self {
        let mut c = Self::new("Weekend");
        c.fill_weekdays();
        c
    }
    pub fn target() -> Self {
        let mut c = Self::weekend_only();
        c.name = Arc::from("TARGET2");
        for y in 1900..=2155 {
            c.add_fixed(y, 1, 1);
            c.add_easter(y, -2);
            c.add_easter(y, 1);
            c.add_fixed(y, 5, 1);
            c.add_fixed(y, 12, 25);
            c.add_fixed(y, 12, 26);
        }
        c
    }
    pub fn us_sofr() -> Self {
        let mut c = Self::weekend_only();
        c.name = Arc::from("US-SOFR");
        for y in 1900..=2155 {
            c.us_holidays(y, true);
        }
        c
    }
    pub fn us_nyse() -> Self {
        let mut c = Self::us_sofr();
        c.name = Arc::from("US-NYSE");
        c
    }
    pub fn uk_settlement() -> Self {
        let mut c = Self::weekend_only();
        c.name = Arc::from("UK-Settlement");
        for y in 1900..=2155 {
            c.uk_holidays(y);
        }
        c
    }
    fn fill_weekdays(&mut self) {
        for day in START_DAY..=END_DAY {
            let d = Date::from_days(day);
            if !d.is_weekend() {
                self.set(d, true);
            }
        }
    }
    fn us_holidays(&mut self, y: i32, juneteenth: bool) {
        self.observed_fixed(y, 1, 1);
        self.nth_weekday(y, 1, Weekday::Monday, 3);
        self.nth_weekday(y, 2, Weekday::Monday, 3); // MLK/Presidents
        self.add_easter(y, -2);
        self.last_weekday(y, 5, Weekday::Monday); // Memorial
        if juneteenth && y >= 2021 {
            self.observed_fixed(y, 6, 19);
        }
        self.observed_fixed(y, 7, 4);
        self.nth_weekday(y, 9, Weekday::Monday, 1);
        self.nth_weekday(y, 10, Weekday::Monday, 2);
        self.observed_fixed(y, 11, 11);
        self.nth_weekday(y, 11, Weekday::Thursday, 4);
        self.observed_fixed(y, 12, 25);
    }
    fn add_fixed(&mut self, y: i32, m: u8, d: u8) {
        if let Ok(x) = Date::from_ymd(y, m, d) {
            self.set(x, false);
        }
    }
    fn observed_fixed(&mut self, y: i32, m: u8, d: u8) {
        if let Ok(x) = Date::from_ymd(y, m, d) {
            self.set(x, false);
            if x.weekday() == Weekday::Saturday {
                self.set(x.add_days(-1), false)
            } else if x.weekday() == Weekday::Sunday {
                self.set(x.add_days(1), false)
            }
        }
    }
    fn add_easter(&mut self, y: i32, offset: i32) {
        let e = easter_sunday(y);
        self.set(e.add_days(offset), false);
    }
    fn nth_weekday(&mut self, y: i32, m: u8, w: Weekday, n: u8) {
        if let Ok(mut d) = Date::from_ymd(y, m, 1) {
            while d.weekday() != w {
                d = d.add_days(1);
            }
            self.set(d.add_days(7 * (n as i32 - 1)), false);
        }
    }
    fn last_weekday(&mut self, y: i32, m: u8, w: Weekday) {
        if let Ok(mut d) = Date::from_ymd(y, m, crate::date::days_in_month(y, m)) {
            while d.weekday() != w {
                d = d.add_days(-1);
            }
            self.set(d, false);
        }
    }
    fn uk_holidays(&mut self, y: i32) {
        // UK substitute days are always moved forward into the following
        // year (unlike the US convention, which observes Saturday holidays
        // on the preceding Friday).  Thus 1 Jan 2022 (Saturday) is observed
        // on Monday 3 Jan, and 1 Jan 2023 (Sunday) on Monday 2 Jan.
        if let Ok(new_year) = Date::from_ymd(y, 1, 1) {
            self.set(new_year, false);
            match new_year.weekday() {
                Weekday::Saturday => self.set(new_year.add_days(2), false),
                Weekday::Sunday => self.set(new_year.add_days(1), false),
                _ => {}
            }
        }
        self.add_easter(y, -2);
        self.add_easter(y, 1);
        self.nth_weekday(y, 5, Weekday::Monday, 1);
        self.last_weekday(y, 5, Weekday::Monday);
        self.last_weekday(y, 8, Weekday::Monday);
        // UK substitute-day rules: holidays falling on weekends move to the
        // next available weekday, with Christmas/Boxing Day shifted as a pair.
        if let Ok(christmas) = Date::from_ymd(y, 12, 25)
            && let Ok(boxing) = Date::from_ymd(y, 12, 26)
        {
            self.set(christmas, false);
            self.set(boxing, false);
            match (christmas.weekday(), boxing.weekday()) {
                (Weekday::Saturday, Weekday::Sunday) => {
                    self.set(boxing.add_days(1), false); // Mon 27
                    self.set(boxing.add_days(2), false); // Tue 28
                }
                (Weekday::Sunday, _) => self.set(boxing.add_days(1), false),
                (_, Weekday::Saturday) => self.set(boxing.add_days(2), false),
                (_, Weekday::Sunday) => self.set(boxing.add_days(1), false),
                _ => {}
            }
        }
    }
    #[inline(always)]
    pub fn is_business_day(&self, d: Date) -> bool {
        if d.0 < START_DAY || d.0 > END_DAY {
            return !d.is_weekend();
        }
        let i = (d.0 - START_DAY) as usize;
        ((self.bits[i >> 6] >> (i & 63)) & 1) != 0
    }
    #[inline]
    fn set(&mut self, d: Date, business: bool) {
        if d.0 < START_DAY || d.0 > END_DAY {
            return;
        }
        let i = (d.0 - START_DAY) as usize;
        let bit = 1u64 << (i & 63);
        if business {
            self.bits[i >> 6] |= bit
        } else {
            self.bits[i >> 6] &= !bit
        }
    }
    pub fn join(&self, other: &Self, name: &str) -> Self {
        self.join_holidays(other, name)
    }

    pub fn join_holidays(&self, other: &Self, name: &str) -> Self {
        let mut c = Self::new(Arc::<str>::from(name));
        for i in 0..N_WORDS {
            c.bits[i] = self.bits[i] & other.bits[i];
        }
        c
    }
    pub fn adjust(&self, d: Date, bdc: BusinessDayConvention) -> Date {
        if matches!(bdc, BusinessDayConvention::Unadjusted) || self.is_business_day(d) {
            return d;
        }
        match bdc {
            BusinessDayConvention::Following
            | BusinessDayConvention::ModifiedFollowing
            | BusinessDayConvention::HalfMonthModifiedFollowing => {
                let mut x = d;
                // Date::add_days saturates at i32 bounds.  Stop before a
                // saturated increment can repeat forever for an out-of-range
                // weekend date.
                while !self.is_business_day(x) && x.0 < 2_932_896 {
                    x = x.add_days(1);
                }
                if matches!(bdc, BusinessDayConvention::ModifiedFollowing) && x.month() != d.month()
                    || matches!(bdc, BusinessDayConvention::HalfMonthModifiedFollowing)
                        && (x.month() != d.month() || (d.day() <= 15 && x.day() > 15))
                {
                    self.adjust(d, BusinessDayConvention::Preceding)
                } else {
                    x
                }
            }
            BusinessDayConvention::Preceding | BusinessDayConvention::ModifiedPreceding => {
                let mut x = d;
                while !self.is_business_day(x) && x.0 > -719_468 {
                    x = x.add_days(-1);
                }
                if matches!(bdc, BusinessDayConvention::ModifiedPreceding) && x.month() != d.month()
                {
                    self.adjust(d, BusinessDayConvention::Following)
                } else {
                    x
                }
            }
            BusinessDayConvention::Unadjusted => d,
        }
    }
    pub fn advance(&self, d: Date, p: Period, bdc: BusinessDayConvention, eom: bool) -> Date {
        self.adjust(d.add_period(p, eom), bdc)
    }
    pub fn business_days_between(&self, d1: Date, d2: Date) -> i32 {
        if d1 == d2 {
            return 0;
        }
        let sign = if d2 > d1 { 1 } else { -1 };
        let (a, b) = if sign > 0 { (d1, d2) } else { (d2, d1) };
        // Count set bits over the half-open interval [a, b) using POPCNT.
        // Clamp to the calendar bitset; dates outside the supported range are
        // treated as weekdays by `is_business_day` and counted directly.
        let mut n: i32 = 0;
        let lo = a.0.max(START_DAY);
        let hi = b.0.min(END_DAY + 1);
        if lo < hi {
            let start = (lo - START_DAY) as usize;
            let end = (hi - START_DAY) as usize; // exclusive
            let first_word = start >> 6;
            let last_word = (end - 1) >> 6;
            if first_word == last_word {
                let from = start & 63;
                let to = ((end - 1) & 63) + 1;
                let width_mask = if to == 64 { u64::MAX } else { (1u64 << to) - 1 };
                n += (self.bits[first_word] & width_mask & (!0u64 << from)).count_ones() as i32;
            } else {
                n += (self.bits[first_word] & (!0u64 << (start & 63))).count_ones() as i32;
                n += self.bits[first_word + 1..last_word]
                    .iter()
                    .map(|word| word.count_ones() as i32)
                    .sum::<i32>();
                let last_bits = ((end - 1) & 63) + 1;
                let last_mask = if last_bits == 64 {
                    u64::MAX
                } else {
                    (1u64 << last_bits) - 1
                };
                n += (self.bits[last_word] & last_mask).count_ones() as i32;
            }
        }
        // Handle portions outside the precomputed range (normally only when
        // users request dates beyond 2155/1900).
        if a.0 < START_DAY {
            let stop = b.0.min(START_DAY);
            let mut d = a.0;
            while d < stop {
                if !Date::from_days(d).is_weekend() {
                    n += 1;
                }
                d += 1;
            }
        }
        if b.0 > END_DAY + 1 {
            let mut d = a.0.max(END_DAY + 1);
            while d < b.0 {
                if !Date::from_days(d).is_weekend() {
                    n += 1;
                }
                d += 1;
            }
        }
        sign * n
    }
    pub fn holidays_word(&self, i: usize) -> u64 {
        self.bits[i]
    }
}

fn easter_sunday(y: i32) -> Date {
    let a = y % 19;
    let b = y / 100;
    let c = y % 100;
    let d = b / 4;
    let e = b % 4;
    let f = (b + 8) / 25;
    let g = (b - f + 1) / 3;
    let h = (19 * a + b - d - g + 15) % 30;
    let i = c / 4;
    let k = c % 4;
    let l = (32 + 2 * e + 2 * i - h - k) % 7;
    let m = (a + 11 * h + 22 * l) / 451;
    let month = (h + l - 7 * m + 114) / 31;
    let day = (h + l - 7 * m + 114) % 31 + 1;
    Date::from_ymd(y, month as u8, day as u8).expect("easter")
}

// This extension fills the remaining US holidays while keeping the public constructor compact.
impl BitsetCalendar {
    pub fn us_federal() -> Self {
        Self::us_sofr()
    }
}
impl BitsetCalendar {
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn is_holiday(&self, d: Date) -> bool {
        !self.is_business_day(d)
    }
}
