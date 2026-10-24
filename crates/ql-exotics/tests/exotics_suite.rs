use ql_ad::Dual;
use ql_exotics::*;
use ql_instruments::OptionType;
use ql_vol::black_scholes_price;

fn spec(bt: BarrierType, call: bool, k: f64) -> BarrierOptionSpec<f64> {
    BarrierOptionSpec {
        option_type: if call {
            OptionType::Call
        } else {
            OptionType::Put
        },
        barrier_type: bt,
        spot: 100.,
        strike: k,
        barrier: 95.,
        rebate: 0.,
        rate: 0.03,
        div_yield: 0.01,
        vol: 0.25,
        expiry: 1.,
        monitoring: MonitoringFrequency::Continuous,
    }
}
#[test]
fn barrier_in_out_parity_all_branches() {
    for &call in &[true, false] {
        for &k in &[90., 95., 110.] {
            let v = black_scholes_price(100., k, 0.25, 1., 0.03, 0.01, call);
            let di = price_barrier_option(&spec(BarrierType::DownIn, call, k));
            let do_ = price_barrier_option(&spec(BarrierType::DownOut, call, k));
            let ui = price_barrier_option(&spec(BarrierType::UpIn, call, k));
            let uo = price_barrier_option(&spec(BarrierType::UpOut, call, k));
            assert!(
                (di + do_ - v).abs() < 1e-10,
                "down {call} {k} {}",
                di + do_ - v
            );
            assert!((ui + uo - v).abs() < 1e-10, "up {call} {k} {}", ui + uo - v);
        }
    }
}
#[test]
fn asian_and_heston_variance() {
    let g = price_geometric_asian_continuous(100., 100., 0.03, 0., 0.2, 1., OptionType::Call);
    let a = price_turnbull_wakeman_arithmetic_asian(
        100.,
        100.,
        0.03,
        0.,
        0.2,
        &[0.25, 0.5, 0.75, 1.],
        1.,
        0,
        0.,
        OptionType::Call,
    );
    assert!(g > 0. && a > g);
    let v = heston_analytical_fair_variance(0.09, 1., 0.04, 1.);
    assert!(v > 0.04 && v < 0.09);
}
#[test]
fn barrier_dual_greeks_live() {
    let s = spec(BarrierType::DownOut, true, 110.);
    let (_, g) = barrier_option_greeks_dual5(&s);
    assert!(g.iter().all(|x| x.is_finite()));
    let d = BarrierOptionSpec {
        spot: Dual::<5>::variable(s.spot, 0),
        strike: Dual::from_value(s.strike),
        barrier: Dual::from_value(s.barrier),
        rebate: Dual::from_value(0.),
        rate: Dual::variable(s.rate, 2),
        div_yield: Dual::variable(s.div_yield, 3),
        vol: Dual::variable(s.vol, 1),
        expiry: Dual::variable(s.expiry, 4),
        option_type: s.option_type,
        barrier_type: s.barrier_type,
        monitoring: s.monitoring,
    };
    assert!(price_barrier_option(&d).d.iter().all(|x| x.is_finite()));
    let (_, greeks) = barrier_option_greeks_dual5(&s);
    for (j, &greek) in greeks.iter().enumerate() {
        let h = 1e-5;
        let mut up = s;
        let mut dn = s;
        match j {
            0 => {
                up.spot += h;
                dn.spot -= h;
            }
            1 => {
                up.vol += h;
                dn.vol -= h;
            }
            2 => {
                up.rate += h;
                dn.rate -= h;
            }
            3 => {
                up.div_yield += h;
                dn.div_yield -= h;
            }
            _ => {
                up.expiry += h;
                dn.expiry -= h;
            }
        }
        let fd = (price_barrier_option(&up) - price_barrier_option(&dn)) / (2. * h);
        let fd_theta = if j == 4 { -fd } else { fd };
        assert!(
            (greek - fd_theta).abs() < 1e-6,
            "j={j} ad={} fd={fd_theta}",
            greek
        );
    }
}

#[test]
fn bgk_monitoring_and_variance_swap() {
    let mut s = spec(BarrierType::DownOut, true, 100.);
    let cont = price_barrier_option(&s);
    let mut previous = cont;
    for n in [4, 12, 52, 252, 100_000] {
        s.monitoring = MonitoringFrequency::Discrete { observations: n };
        let p = price_barrier_option(&s);
        assert!(p.is_finite());
        previous = p;
    }
    assert!(previous > 0.0 && cont > 0.0);
    let mut q = Vec::new();
    for i in 1..=101 {
        let k = 50.0 + i as f64;
        q.push(VarianceSwapStrikeQuote {
            strike: k,
            implied_vol: 0.2,
        });
    }
    let (_, vol) = replicate_variance_swap_fair_strike(100., 0.03, 0., 1., 100., &q);
    assert!((vol - 0.2).abs() < 5e-4);
    assert!(volatility_swap_convexity_adjustment(0.04, 0.0001) < 0.2);
}

#[test]
fn phase11_variance_and_corridor_replication() {
    let quotes = (0..21)
        .map(|i| VarianceSwapStrikeQuote {
            strike: 60.0 + i as f64 * 4.0,
            implied_vol: 0.2,
        })
        .collect::<Vec<_>>();
    let base = replicate_variance_swap_fair_strike(100.0, 0.02, 0.0, 1.0, 100.0, &quotes);
    let wide = replicate_corridor_variance_swap_fair_strike(
        100.0, 0.02, 0.0, 1.0, 100.0, 1e-4, 1e4, &quotes,
    );
    assert!((base.0 - wide.0).abs() < 1e-10);
    let small = heston_integrated_variance_var(0.04, 1e-8, 0.04, 0.3, 1.0);
    assert!(small.is_finite() && small >= 0.0);
    let left = replicate_corridor_variance_swap_fair_strike(
        100.0, 0.02, 0.0, 1.0, 100.0, 70.0, 100.0, &quotes,
    )
    .0;
    let right = replicate_corridor_variance_swap_fair_strike(
        100.0, 0.02, 0.0, 1.0, 100.0, 100.0, 140.0, &quotes,
    )
    .0;
    let all = replicate_corridor_variance_swap_fair_strike(
        100.0, 0.02, 0.0, 1.0, 100.0, 70.0, 140.0, &quotes,
    )
    .0;
    assert!((left + right - all).abs() < 1e-8);
}

#[test]
fn phase12_round4_autocallable_price_greeks_and_life_invariants() {
    let s = WorstOfAutocallableSpec {
        notional: 1.0,
        spots: vec![100.25, 100.0],
        vols: vec![0.2, 0.2],
        divs: vec![0.0, 0.0],
        rate: 0.02,
        correlation: vec![vec![1.0, 0.2], vec![0.2, 1.0]],
        observation_times: vec![0.5, 1.0],
        autocall_barriers: vec![1.1, 1.05],
        redemption_bonuses: vec![0.0, 0.0],
        coupon_barriers: vec![0.8, 0.8],
        coupon_rates: vec![0.01, 0.01],
        memory_coupon: true,
        knock_in_barrier: 0.6,
        put_strike: 1.0,
        observe_ki_all_dates: true,
        smoothing_half_width: 0.02,
        smooth_min_temp: 0.01,
    };
    let r = price_worst_of_autocallable_aad(&s, 64, 4, 7).unwrap();
    assert!(r.price.is_finite() && r.stderr.is_finite());
    assert_eq!(r.autocall_probabilities.len(), s.observation_times.len());
    assert!(
        r.autocall_probabilities
            .iter()
            .all(|p| *p >= 0.0 && *p <= 1.0)
    );
    assert!(r.autocall_probabilities.iter().sum::<f64>() <= 1.0 + 1e-12);
    assert!(r.expected_life >= 0.0 && r.expected_life <= 1.0 + 1e-12);
    let h = 1e-4;
    let mut up = s.clone();
    let mut dn = s.clone();
    up.spots[0] += h;
    dn.spots[0] -= h;
    let fd_spot = (price_worst_of_autocallable(&up, 64, 7).unwrap()
        - price_worst_of_autocallable(&dn, 64, 7).unwrap())
        / (2.0 * h);
    assert!((r.spot_deltas[0] - fd_spot).abs() < 1e-6);
    up = s.clone();
    dn = s.clone();
    up.vols[0] += h;
    dn.vols[0] -= h;
    let fd_vega = (price_worst_of_autocallable(&up, 64, 7).unwrap()
        - price_worst_of_autocallable(&dn, 64, 7).unwrap())
        / (2.0 * h);
    assert!((r.vegas[0] - fd_vega).abs() < 1e-5);
    up = s.clone();
    dn = s.clone();
    up.correlation[0][1] += h;
    up.correlation[1][0] += h;
    dn.correlation[0][1] -= h;
    dn.correlation[1][0] -= h;
    let fd_cega = (price_worst_of_autocallable(&up, 64, 7).unwrap()
        - price_worst_of_autocallable(&dn, 64, 7).unwrap())
        / (2.0 * h);
    assert!((r.corr_cega[1][0] - fd_cega).abs() < 1e-5);
}
