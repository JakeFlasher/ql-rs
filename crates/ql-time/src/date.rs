use serde::{Deserialize, Serialize};
use thiserror::Error;

#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub struct Date(pub i32);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Error)]
pub enum TimeError {
    #[error("invalid date: {0:04}-{1:02}-{2:02}")]
    InvalidDate(i32, u8, u8),
    #[error("invalid month {0}")]
    InvalidMonth(u8),
    #[error("date out of supported calendar range")]
    OutOfRange,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Weekday {
    Monday,
    Tuesday,
    Wednesday,
    Thursday,
    Friday,
    Saturday,
    Sunday,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Frequency {
    Once,
    NoFrequency,
    Annual,
    Semiannual,
    Biannual,
    Quarterly,
    EveryFourthMonth,
    Bimonthly,
    Monthly,
    Weekly,
    Daily,
}

impl Frequency {
    pub fn months(self) -> i32 {
        match self {
            Self::Once | Self::NoFrequency => 0,
            Self::Annual => 12,
            Self::Semiannual | Self::Biannual => 6,
            Self::Quarterly => 3,
            Self::EveryFourthMonth => 4,
            Self::Bimonthly => 2,
            Self::Monthly => 1,
            Self::Weekly => 0,
            Self::Daily => 0,
        }
    }
    pub fn periods_per_year(self) -> i32 {
        match self {
            Self::Once | Self::NoFrequency => 1,
            Self::Annual => 1,
            Self::Semiannual | Self::Biannual => 2,
            Self::Quarterly => 4,
            Self::EveryFourthMonth => 3,
            Self::Bimonthly => 6,
            Self::Monthly => 12,
            Self::Weekly => 52,
            Self::Daily => 365,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Period {
    Days(i32),
    Weeks(i32),
    Months(i32),
    Years(i32),
    Frequency(i32, Frequency),
}
impl Period {
    pub fn months(self) -> i32 {
        match self {
            Self::Months(x) => x,
            Self::Years(x) => x.saturating_mul(12),
            Self::Frequency(n, f) => n.saturating_mul(f.months()),
            _ => 0,
        }
    }
}

impl Date {
    pub const EPOCH: Self = Self(0);
    pub const MIN: Self = Self(-719_468); // 0000-03-01-ish, practical lower bound
    pub const MAX: Self = Self(i32::MAX);
    pub const fn from_days(days: i32) -> Self {
        Self(days)
    }
    pub const fn days(self) -> i32 {
        self.0
    }
    pub const fn to_epoch_days(self) -> i32 {
        self.0
    }
    pub fn from_ymd(year: i32, month: u8, day: u8) -> Result<Self, TimeError> {
        if !(1..=12).contains(&month) {
            return Err(TimeError::InvalidMonth(month));
        }
        let dim = days_in_month(year, month);
        if day == 0 || day > dim {
            return Err(TimeError::InvalidDate(year, month, day));
        }
        let days = days_from_civil(year, month as i32, day as i32);
        if days < i32::MIN as i64 || days > i32::MAX as i64 {
            return Err(TimeError::OutOfRange);
        }
        Ok(Self(days as i32))
    }
    pub fn ymd(self) -> (i32, u8, u8) {
        civil_from_days(self.0)
    }
    pub fn year_month_day(self) -> (i32, u8, u8) {
        self.ymd()
    }
    pub fn year(self) -> i32 {
        self.ymd().0
    }
    pub fn month(self) -> u8 {
        self.ymd().1
    }
    pub fn day(self) -> u8 {
        self.ymd().2
    }
    pub fn weekday(self) -> Weekday {
        match self.0.rem_euclid(7) {
            0 => Weekday::Thursday,
            1 => Weekday::Friday,
            2 => Weekday::Saturday,
            3 => Weekday::Sunday,
            4 => Weekday::Monday,
            5 => Weekday::Tuesday,
            _ => Weekday::Wednesday,
        }
    }
    pub fn is_weekend(self) -> bool {
        matches!(self.weekday(), Weekday::Saturday | Weekday::Sunday)
    }
    pub const fn is_leap_year(year: i32) -> bool {
        year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)
    }
    pub fn end_of_month(self) -> Self {
        Self::from_ymd(
            self.year(),
            self.month(),
            days_in_month(self.year(), self.month()),
        )
        .expect("valid")
    }
    pub const fn add_days(self, days: i32) -> Self {
        Self(self.0.saturating_add(days))
    }
    pub fn add_months(self, months: i32, eom: bool) -> Self {
        let (y, m, d) = self.ymd();
        let z = (y as i64) * 12 + (m as i64 - 1) + months as i64;
        let ny64 = z.div_euclid(12);
        if ny64 < i32::MIN as i64 || ny64 > i32::MAX as i64 {
            return if z < 0 { Self::MIN } else { Self::MAX };
        }
        let ny = ny64 as i32;
        let nm = (z.rem_euclid(12) + 1) as u8;
        let nd = if eom && d == days_in_month(y, m) {
            days_in_month(ny, nm)
        } else {
            d.min(days_in_month(ny, nm))
        };
        Self::from_ymd(ny, nm, nd).unwrap_or(if z < 0 { Self::MIN } else { Self::MAX })
    }
    pub fn imm_date(year: i32, month: u8) -> Result<Self, TimeError> {
        crate::imm_date(year, month)
    }
    pub fn cds_date(year: i32, month: u8) -> Result<Self, TimeError> {
        crate::cds_date(year, month)
    }
    pub fn add_period(self, p: Period, eom: bool) -> Self {
        match p {
            Period::Days(n) => self.add_days(n),
            Period::Weeks(n) => self.add_days(n.saturating_mul(7)),
            Period::Months(n) => self.add_months(n, eom),
            Period::Years(n) => self.add_months(n.saturating_mul(12), eom),
            Period::Frequency(n, f) => {
                if f.months() != 0 {
                    self.add_months(n.saturating_mul(f.months()), eom)
                } else {
                    self.add_days(match f {
                        Frequency::Weekly => n.saturating_mul(7),
                        Frequency::Daily => n,
                        _ => n,
                    })
                }
            }
        }
    }
}

pub fn days_in_month(year: i32, month: u8) -> u8 {
    match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 => {
            if Date::is_leap_year(year) {
                29
            } else {
                28
            }
        }
        _ => 0,
    }
}
pub fn is_leap_year(year: i32) -> bool {
    Date::is_leap_year(year)
}
fn days_from_civil(y: i32, m: i32, d: i32) -> i64 {
    // Keep every intermediate in i64.  A valid i32 year can be far outside
    // Date's representable epoch-day range, and doing the arithmetic in i32
    // used to wrap before the range check could run.
    let y = y as i64 - if m <= 2 { 1 } else { 0 };
    let era = (if y >= 0 { y } else { y - 399 }).div_euclid(400);
    let yoe = y - era * 400;
    let mp = m as i64 + if m > 2 { -3 } else { 9 };
    let doy = (153 * mp + 2).div_euclid(5) + d as i64 - 1;
    let doe = yoe * 365 + yoe.div_euclid(4) - yoe.div_euclid(100) + doy;
    era * 146097 + doe - 719468
}
fn civil_from_days(z: i32) -> (i32, u8, u8) {
    let z = z as i64 + 719468;
    let era = (if z >= 0 { z } else { z - 146096 }).div_euclid(146097);
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096).div_euclid(365);
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2).div_euclid(153);
    let d = doy - (153 * mp + 2).div_euclid(5) + 1;
    let m = mp + if mp < 10 { 3 } else { -9 };
    ((y + if m <= 2 { 1 } else { 0 }) as i32, m as u8, d as u8)
}

pub fn imm_date(year: i32, month: u8) -> Result<Date, TimeError> {
    let first = Date::from_ymd(year, month, 1)?;
    let wd: i32 = match first.weekday() {
        Weekday::Monday => 0,
        Weekday::Tuesday => 1,
        Weekday::Wednesday => 2,
        Weekday::Thursday => 3,
        Weekday::Friday => 4,
        Weekday::Saturday => 5,
        Weekday::Sunday => 6,
    };
    let day = 1 + (2 - wd).rem_euclid(7) + 14;
    Date::from_ymd(year, month, day as u8)
}
pub fn cds_date(year: i32, month: u8) -> Result<Date, TimeError> {
    if matches!(month, 3 | 6 | 9 | 12) {
        Date::from_ymd(year, month, 20)
    } else {
        Err(TimeError::InvalidMonth(month))
    }
}
pub fn cds_imm_date(year: i32, month: u8) -> Result<Date, TimeError> {
    cds_date(year, month)
}
impl std::fmt::Display for Date {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (y, m, d) = self.ymd();
        write!(f, "{y:04}-{m:02}-{d:02}")
    }
}
impl From<i32> for Date {
    fn from(v: i32) -> Self {
        Self(v)
    }
}
impl From<Date> for i32 {
    fn from(v: Date) -> Self {
        v.0
    }
}
