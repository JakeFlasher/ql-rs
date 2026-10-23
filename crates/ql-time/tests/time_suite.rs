use ql_time::*;
#[test]
fn civil_round_trip_and_imm() {
    for y in 1900..=2155 {
        for m in 1..=12 {
            for d in 1..=days_in_month(y, m) {
                let x = Date::from_ymd(y, m, d).unwrap();
                assert_eq!(x.ymd(), (y, m, d));
            }
        }
    }
    assert_eq!(imm_date(2024, 3).unwrap().ymd(), (2024, 3, 20));
    assert_eq!(cds_date(2024, 3).unwrap().ymd(), (2024, 3, 20));
}
#[test]
fn holiday_rules() {
    let uk = BitsetCalendar::uk_settlement();
    assert!(!uk.is_business_day(Date::from_ymd(2021, 12, 27).unwrap()));
    assert!(!uk.is_business_day(Date::from_ymd(2021, 12, 28).unwrap()));
    assert!(uk.is_business_day(Date::from_ymd(2021, 12, 24).unwrap()));
    let ny = BitsetCalendar::us_nyse();
    assert!(!ny.is_business_day(Date::from_ymd(2024, 3, 29).unwrap()));
    assert!(!ny.is_business_day(Date::from_ymd(2024, 3, 30).unwrap()));
}
#[test]
fn popcnt_and_daycounts() {
    let c = BitsetCalendar::weekend_only();
    assert_eq!(
        c.business_days_between(
            Date::from_ymd(2024, 1, 1).unwrap(),
            Date::from_ymd(2024, 1, 8).unwrap()
        ),
        5
    );
    assert!(
        (DayCounter::Act360.year_fraction(
            Date::from_ymd(2024, 1, 1).unwrap(),
            Date::from_ymd(2024, 7, 1).unwrap(),
            None,
            None
        ) - 182.0 / 360.0)
            .abs()
            < 1e-14
    );
}
#[test]
fn schedule_anchor() {
    let cal = std::sync::Arc::new(BitsetCalendar::weekend_only());
    let s = ScheduleBuilder::new(
        Date::from_ymd(2021, 1, 31).unwrap(),
        Date::from_ymd(2022, 7, 31).unwrap(),
        Period::Months(1),
        cal,
    )
    .rule(DateGenerationRule::Forward)
    .convention(BusinessDayConvention::Unadjusted)
    .end_of_month(true)
    .build();
    assert_eq!(s.dates()[1].day(), 28);
    assert_eq!(s.dates()[2].day(), 31);
    assert_eq!(s.dates()[2].day(), 31);
}
