use ql_rates::{CmsSpreadSpec, LognormalCmsSpreadPricer};

fn spec() -> CmsSpreadSpec<f64> {
    CmsSpreadSpec {
        fixing_time: 1.0,
        payment_time: 1.0,
        accrual_factor: 0.5,
        discount_df: 0.98,
        gearing1: 1.0,
        gearing2: 1.0,
        spread: 0.0,
        swap_rate1: 0.04,
        swap_rate2: 0.02,
        vol1: 0.2,
        vol2: 0.25,
        rho: 0.2,
        shift1: 0.0,
        shift2: 0.0,
        mean_reversion: 0.0,
        annuity1: 1.0,
        annuity2: 1.0,
        swap_tenor1: 10.0,
        swap_tenor2: 5.0,
    }
}

#[test]
fn cms_spread_parity_and_collar() {
    let s = spec();
    let k = 0.01;
    let cap = LognormalCmsSpreadPricer::caplet_pv(&s, k);
    let floor = LognormalCmsSpreadPricer::floorlet_pv(&s, k);
    let (f, pv) = LognormalCmsSpreadPricer::coupon_rate_and_pv(&s);
    assert!((cap - floor - 0.98 * 0.5 * (f - k)).abs() < 1e-12);
    let collar = LognormalCmsSpreadPricer::collar_pv(&s, 0.03, 0.0);
    assert!(
        (collar
            - (LognormalCmsSpreadPricer::caplet_pv(&s, 0.03)
                - LognormalCmsSpreadPricer::floorlet_pv(&s, 0.0)))
        .abs()
            < 1e-12
    );
    assert!(pv > 0.0);
}
