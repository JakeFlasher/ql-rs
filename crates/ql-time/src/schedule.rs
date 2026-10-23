use crate::{BitsetCalendar, BusinessDayConvention, Date, Frequency, Period};
use smallvec::SmallVec;
use std::sync::Arc;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DateGenerationRule {
    Backward,
    Forward,
    Zero,
    ThirdWednesday,
    Twentieth,
    TwentiethImm,
}
#[derive(Debug, Clone)]
pub struct Schedule {
    pub dates: SmallVec<[Date; 64]>,
    pub is_regular: SmallVec<[bool; 64]>,
    pub convention: BusinessDayConvention,
    pub calendar: Arc<BitsetCalendar>,
    pub end_of_month: bool,
}
impl Schedule {
    pub fn dates(&self) -> &[Date] {
        &self.dates
    }
    pub fn is_regular(&self) -> &[bool] {
        &self.is_regular
    }
    pub fn len(&self) -> usize {
        self.dates.len()
    }
    pub fn is_empty(&self) -> bool {
        self.dates.is_empty()
    }
}
#[derive(Debug, Clone)]
pub struct ScheduleBuilder {
    effective_date: Date,
    termination_date: Date,
    tenor: Period,
    calendar: Arc<BitsetCalendar>,
    convention: BusinessDayConvention,
    termination_convention: Option<BusinessDayConvention>,
    rule: DateGenerationRule,
    end_of_month: bool,
    first_date: Option<Date>,
    next_to_last_date: Option<Date>,
}
impl ScheduleBuilder {
    pub fn new(
        effective_date: Date,
        termination_date: Date,
        tenor: Period,
        calendar: Arc<BitsetCalendar>,
    ) -> Self {
        Self {
            effective_date,
            termination_date,
            tenor,
            calendar,
            convention: BusinessDayConvention::Following,
            termination_convention: None,
            rule: DateGenerationRule::Backward,
            end_of_month: false,
            first_date: None,
            next_to_last_date: None,
        }
    }
    pub fn convention(mut self, x: BusinessDayConvention) -> Self {
        self.convention = x;
        self
    }
    pub fn termination_convention(mut self, x: BusinessDayConvention) -> Self {
        self.termination_convention = Some(x);
        self
    }
    pub fn rule(mut self, x: DateGenerationRule) -> Self {
        self.rule = x;
        self
    }
    pub fn date_generation_rule(self, x: DateGenerationRule) -> Self {
        self.rule(x)
    }
    pub fn end_of_month(mut self, x: bool) -> Self {
        self.end_of_month = x;
        self
    }
    pub fn first_date(mut self, x: Option<Date>) -> Self {
        self.first_date = x;
        self
    }
    pub fn next_to_last_date(mut self, x: Option<Date>) -> Self {
        self.next_to_last_date = x;
        self
    }
    pub fn build(self) -> Schedule {
        assert!(
            !matches!(
                self.tenor,
                Period::Days(0)
                    | Period::Weeks(0)
                    | Period::Months(0)
                    | Period::Years(0)
                    | Period::Frequency(0, _)
            ),
            "schedule tenor must be non-zero"
        );
        assert!(
            self.effective_date <= self.termination_date,
            "effective date must not be after termination date"
        );
        let mut raw: SmallVec<[Date; 64]> = SmallVec::new();
        match self.rule {
            DateGenerationRule::Zero => {
                raw.push(self.effective_date);
                raw.push(self.termination_date);
            }
            DateGenerationRule::Forward
            | DateGenerationRule::ThirdWednesday
            | DateGenerationRule::Twentieth
            | DateGenerationRule::TwentiethImm => {
                raw.push(self.effective_date);
                if let Some(x) = self.first_date
                    && x > self.effective_date
                {
                    raw.push(x);
                }
                let mut k = 1;
                loop {
                    let n = self
                        .effective_date
                        .add_period(scale_period(self.tenor, k), self.end_of_month);
                    if n >= self.termination_date {
                        break;
                    }
                    // A zero (or otherwise non-advancing) tenor cannot
                    // generate a schedule.  Without this guard the loop
                    // would run forever because `scale_period` keeps
                    // producing the effective date (e.g. `Period::Days(0)`
                    // or `Frequency::NoFrequency`).
                    if n <= raw.last().copied().unwrap_or(self.effective_date) {
                        break;
                    }
                    if raw.last().copied() != Some(n) {
                        raw.push(n);
                    }
                    k += 1;
                }
                if let Some(x) = self.next_to_last_date
                    && x > raw.last().copied().unwrap_or(self.effective_date)
                    && x < self.termination_date
                {
                    raw.push(x);
                }
                raw.push(self.termination_date);
            }
            DateGenerationRule::Backward => {
                let mut rev: SmallVec<[Date; 64]> = SmallVec::new();
                rev.push(self.termination_date);
                if let Some(x) = self.next_to_last_date
                    && x < self.termination_date
                {
                    rev.push(x);
                }
                let mut k = 1;
                loop {
                    let p = self
                        .termination_date
                        .add_period(scale_period(self.tenor, -k), self.end_of_month);
                    if p <= self.effective_date {
                        break;
                    }
                    // Guard against a non-advancing tenor for the backward
                    // generation rule as well; otherwise a zero tenor makes
                    // this loop non-terminating.
                    if p >= rev.last().copied().unwrap_or(self.termination_date) {
                        break;
                    }
                    if rev.last().copied() != Some(p) {
                        rev.push(p);
                    }
                    k += 1;
                }
                if let Some(x) = self.first_date
                    && x > self.effective_date
                    && x < self.termination_date
                {
                    rev.push(x);
                }
                rev.push(self.effective_date);
                rev.reverse();
                raw = rev;
            }
        }
        // Explicit stub dates can be supplied together with a backward
        // generation rule.  Normalize them before adjustment so the returned
        // schedule is strictly chronological and contains no duplicate node.
        raw.sort_unstable();
        raw.dedup();
        let mut dates: SmallVec<[Date; 64]> = SmallVec::new();
        for (i, d) in raw.iter().enumerate() {
            let conv = if i + 1 == raw.len() {
                self.termination_convention.unwrap_or(self.convention)
            } else {
                self.convention
            };
            let mut x = *d;
            if matches!(self.rule, DateGenerationRule::ThirdWednesday) {
                if let Ok(im) = crate::imm_date(d.year(), d.month()) {
                    x = im;
                }
            } else if matches!(self.rule, DateGenerationRule::Twentieth) {
                if let Ok(z) = Date::from_ymd(d.year(), d.month(), 20) {
                    x = z;
                }
            } else if matches!(self.rule, DateGenerationRule::TwentiethImm)
                && let Ok(z) = crate::cds_date(d.year(), d.month())
            {
                x = z;
            }
            dates.push(self.calendar.adjust(x, conv));
        }
        let mut regular = SmallVec::new();
        for i in 0..raw.len().saturating_sub(1) {
            regular.push(raw[i].add_period(self.tenor, self.end_of_month) == raw[i + 1]);
        }
        Schedule {
            dates,
            is_regular: regular,
            convention: self.convention,
            calendar: self.calendar,
            end_of_month: self.end_of_month,
        }
    }
}
fn scale_period(p: Period, n: i32) -> Period {
    match p {
        Period::Days(v) => Period::Days(v * n),
        Period::Weeks(v) => Period::Weeks(v * n),
        Period::Months(v) => Period::Months(v * n),
        Period::Years(v) => Period::Years(v * n),
        Period::Frequency(v, f) => Period::Frequency(v * n, f),
    }
}
impl From<Frequency> for Period {
    fn from(f: Frequency) -> Self {
        match f {
            Frequency::Annual => Self::Years(1),
            Frequency::Semiannual | Frequency::Biannual => Self::Months(6),
            Frequency::Quarterly => Self::Months(3),
            Frequency::EveryFourthMonth => Self::Months(4),
            Frequency::Bimonthly => Self::Months(2),
            Frequency::Monthly => Self::Months(1),
            Frequency::Weekly => Self::Weeks(1),
            Frequency::Daily => Self::Days(1),
            Frequency::Once | Frequency::NoFrequency => Self::Years(1),
        }
    }
}
