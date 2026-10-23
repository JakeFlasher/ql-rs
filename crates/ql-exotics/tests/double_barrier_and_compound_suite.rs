use ql_ad::Dual;
use ql_exotics::*;
use ql_instruments::OptionType;
use ql_vol::black_scholes_price;

fn double_spec(
    option_type: OptionType,
    barrier_lo: f64,
    barrier_hi: f64,
    barrier_type: DoubleBarrierType,
) -> DoubleBarrierOptionSpec<f64> {
    DoubleBarrierOptionSpec {
        option_type,
        barrier_type,
        spot: 100.0,
        strike: 100.0,
        barrier_lo,
        barrier_hi,
        rate: 0.03,
        div_yield: 0.01,
        vol: 0.25,
        expiry: 1.0,
        series_terms: 20,
    }
}

fn compound_spec(mother_type: OptionType, daughter_type: OptionType) -> CompoundOptionSpec<f64> {
    CompoundOptionSpec {
        mother_type,
        daughter_type,
        spot: 100.0,
        strike_mother: 5.0,
        strike_daughter: 100.0,
        rate: 0.03,
        div_yield: 0.01,
        vol: 0.20,
        expiry_mother: 0.50,
        expiry_daughter: 1.0,
    }
}

#[test]
fn double_barrier_in_out_parity_and_single_barrier_limits() {
    for &option_type in &[OptionType::Call, OptionType::Put] {
        let vanilla =
            black_scholes_price(100.0, 100.0, 0.25, 1.0, 0.03, 0.01, option_type.is_call());
        let ko = price_double_barrier_option(&double_spec(
            option_type,
            90.0,
            120.0,
            DoubleBarrierType::KnockOut,
        ));
        let ki = price_double_barrier_option(&double_spec(
            option_type,
            90.0,
            120.0,
            DoubleBarrierType::KnockIn,
        ));
        assert!((ko + ki - vanilla).abs() < 1e-12);

        let high = price_double_barrier_option(&double_spec(
            option_type,
            90.0,
            1000.0,
            DoubleBarrierType::KnockOut,
        ));
        let down = price_barrier_option(&BarrierOptionSpec {
            option_type,
            barrier_type: BarrierType::DownOut,
            spot: 100.0,
            strike: 100.0,
            barrier: 90.0,
            rebate: 0.0,
            rate: 0.03,
            div_yield: 0.01,
            vol: 0.25,
            expiry: 1.0,
            monitoring: MonitoringFrequency::Continuous,
        });
        assert!((high - down).abs() < 1e-10);

        let low = price_double_barrier_option(&double_spec(
            option_type,
            0.01,
            115.0,
            DoubleBarrierType::KnockOut,
        ));
        let up = price_barrier_option(&BarrierOptionSpec {
            option_type,
            barrier_type: BarrierType::UpOut,
            spot: 100.0,
            strike: 100.0,
            barrier: 115.0,
            rebate: 0.0,
            rate: 0.03,
            div_yield: 0.01,
            vol: 0.25,
            expiry: 1.0,
            monitoring: MonitoringFrequency::Continuous,
        });
        assert!((low - up).abs() < 1e-10);
    }
}

#[test]
fn double_barrier_dual_greeks_match_central_fd() {
    let base = double_spec(OptionType::Call, 85.0, 120.0, DoubleBarrierType::KnockOut);
    let dual = DoubleBarrierOptionSpec {
        option_type: base.option_type,
        barrier_type: base.barrier_type,
        spot: Dual::<5>::variable(base.spot, 0),
        strike: Dual::from_value(base.strike),
        barrier_lo: Dual::from_value(base.barrier_lo),
        barrier_hi: Dual::from_value(base.barrier_hi),
        rate: Dual::<5>::variable(base.rate, 2),
        div_yield: Dual::<5>::variable(base.div_yield, 3),
        vol: Dual::<5>::variable(base.vol, 1),
        expiry: Dual::<5>::variable(base.expiry, 4),
        series_terms: base.series_terms,
    };
    let value = price_double_barrier_option(&dual);
    assert!(value.val.is_finite() && value.d.iter().all(|x| x.is_finite()));
    for (lane, &ad_lane) in value.d.iter().enumerate() {
        let h = 1e-5;
        let mut up = base;
        let mut down = base;
        match lane {
            0 => {
                up.spot += h;
                down.spot -= h;
            }
            1 => {
                up.vol += h;
                down.vol -= h;
            }
            2 => {
                up.rate += h;
                down.rate -= h;
            }
            3 => {
                up.div_yield += h;
                down.div_yield -= h;
            }
            _ => {
                up.expiry += h;
                down.expiry -= h;
            }
        }
        let fd =
            (price_double_barrier_option(&up) - price_double_barrier_option(&down)) / (2.0 * h);
        assert!(
            (ad_lane - fd).abs() < 1e-6,
            "lane {lane}: {ad_lane} vs {fd}"
        );
    }
}

#[test]
fn compound_put_call_parity_and_complex_chooser_reduction() {
    for &daughter_type in &[OptionType::Call, OptionType::Put] {
        let call = price_compound_option(&compound_spec(OptionType::Call, daughter_type));
        let put = price_compound_option(&compound_spec(OptionType::Put, daughter_type));
        let daughter =
            black_scholes_price(100.0, 100.0, 0.20, 1.0, 0.03, 0.01, daughter_type.is_call());
        let rhs = daughter - 5.0 * (-0.03_f64 * 0.5).exp();
        assert!((call - put - rhs).abs() < 1e-10);
    }

    let complex = ComplexChooserSpec {
        spot: 100.0,
        strike_call: 100.0,
        strike_put: 100.0,
        choosing_time: 0.5,
        call_expiry: 1.0,
        put_expiry: 1.0,
        rate: 0.03,
        div_yield: 0.01,
        vol: 0.20,
    };
    let expected = price_simple_chooser(100.0, 100.0, 0.03, 0.01, 0.20, 1.0, 0.5);
    assert!((price_complex_chooser(&complex) - expected).abs() < 1e-10);
}

#[test]
fn compound_and_complex_chooser_dual_greeks_match_fd() {
    let base = compound_spec(OptionType::Call, OptionType::Call);
    let dual = CompoundOptionSpec {
        mother_type: base.mother_type,
        daughter_type: base.daughter_type,
        spot: Dual::<5>::variable(base.spot, 0),
        strike_mother: Dual::<5>::variable(base.strike_mother, 1),
        strike_daughter: Dual::from_value(base.strike_daughter),
        rate: Dual::<5>::variable(base.rate, 2),
        div_yield: Dual::from_value(base.div_yield),
        vol: Dual::<5>::variable(base.vol, 3),
        expiry_mother: Dual::<5>::variable(base.expiry_mother, 4),
        expiry_daughter: Dual::from_value(base.expiry_daughter),
    };
    let ad = price_compound_option(&dual);
    for (lane, &ad_lane) in ad.d.iter().enumerate() {
        let h = 1e-5;
        let mut up = base;
        let mut down = base;
        match lane {
            0 => {
                up.spot += h;
                down.spot -= h;
            }
            1 => {
                up.strike_mother += h;
                down.strike_mother -= h;
            }
            2 => {
                up.rate += h;
                down.rate -= h;
            }
            3 => {
                up.vol += h;
                down.vol -= h;
            }
            _ => {
                up.expiry_mother += h;
                down.expiry_mother -= h;
            }
        }
        let fd = (price_compound_option(&up) - price_compound_option(&down)) / (2.0 * h);
        assert!((ad_lane - fd).abs() < 1e-6, "compound lane {lane}");
    }

    let complex = ComplexChooserSpec {
        spot: 100.0,
        strike_call: 100.0,
        strike_put: 100.0,
        choosing_time: 0.5,
        call_expiry: 1.0,
        put_expiry: 1.0,
        rate: 0.03,
        div_yield: 0.01,
        vol: 0.20,
    };
    let c = ComplexChooserSpec {
        spot: Dual::<4>::variable(100.0, 0),
        strike_call: Dual::from_value(100.0),
        strike_put: Dual::from_value(100.0),
        choosing_time: Dual::from_value(0.5),
        call_expiry: Dual::from_value(1.0),
        put_expiry: Dual::from_value(1.0),
        rate: Dual::variable(0.03, 1),
        div_yield: Dual::variable(0.01, 2),
        vol: Dual::variable(0.20, 3),
    };
    let ad_c = price_complex_chooser(&c);
    for (lane, &ad_lane) in ad_c.d.iter().enumerate() {
        let h = 1e-5;
        let mut up = complex;
        let mut down = complex;
        match lane {
            0 => {
                up.spot += h;
                down.spot -= h;
            }
            1 => {
                up.rate += h;
                down.rate -= h;
            }
            2 => {
                up.div_yield += h;
                down.div_yield -= h;
            }
            _ => {
                up.vol += h;
                down.vol -= h;
            }
        }
        let fd = (price_complex_chooser(&up) - price_complex_chooser(&down)) / (2.0 * h);
        assert!((ad_lane - fd).abs() < 1e-6, "chooser lane {lane}");
    }
}
