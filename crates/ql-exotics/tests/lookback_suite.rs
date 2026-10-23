use ql_ad::{Dual, Real};
use ql_exotics::*;
use ql_instruments::OptionType;
use ql_vol::black_scholes_price;

#[test]
fn lookback_lambda_limit_and_vanilla_bound() {
    for rate in [0.03, 0.0] {
        let floating =
            price_floating_lookback_continuous(100.0, 90.0, rate, 0.0, 0.2, 1.0, OptionType::Call);
        let fixed = price_fixed_lookback_continuous(
            100.0,
            90.0,
            100.0,
            rate,
            0.0,
            0.2,
            1.0,
            OptionType::Call,
        );
        let vanilla = black_scholes_price(100.0, 100.0, 0.2, 1.0, rate, 0.0, true);
        assert!(floating > vanilla);
        assert!(fixed > vanilla);
        assert!(floating.is_finite() && fixed.is_finite());
    }
}

#[test]
fn quantlib_reference_fixed_value() {
    let x = price_fixed_lookback_continuous(
        120.0,
        100.0,
        0.0,
        0.10,
        0.06,
        0.30,
        0.50,
        OptionType::Call,
    );
    assert!(x.is_finite());
    for &(m, k, call, expected) in &[
        (110.0, 95.0, true, 31.163_867_260_487_72),
        (100.0, 105.0, true, 20.604_933_499_041_01),
        (90.0, 105.0, false, 26.573_222_724_055_15),
        (100.0, 95.0, false, 16.226_584_162_874_66),
    ] {
        let v = price_fixed_lookback_continuous(
            100.0,
            m,
            k,
            0.05,
            0.05,
            0.3,
            1.0,
            if call {
                OptionType::Call
            } else {
                OptionType::Put
            },
        );
        assert!((v - expected).abs() < 1e-10, "{m} {k} {call}: {v}");
    }
}

#[test]
fn chooser_and_quanto_are_ad_safe() {
    let c = price_simple_chooser(100.0, 100.0, 0.03, 0.01, 0.2, 1.0, 1e-8);
    let call = black_scholes_price(100.0, 100.0, 0.2, 1.0, 0.03, 0.01, true);
    let put = black_scholes_price(100.0, 100.0, 0.2, 1.0, 0.03, 0.01, false);
    assert!((c - call.max(put)).abs() < 1e-8);
    let x = price_quanto_option(
        Dual::<5>::variable(100.0, 0),
        Dual::from_f64(100.0),
        Dual::<5>::variable(0.2, 1),
        Dual::variable(1.0, 2),
        Dual::variable(0.03, 3),
        Dual::from_f64(0.01),
        Dual::variable(0.15, 4),
        Dual::from_f64(-0.25),
        true,
    );
    assert!(x.val().is_finite());
    assert!(x.d.iter().all(|v| v.is_finite()));
}

#[test]
fn bgk_discrete_lookback_converges_monotonically() {
    let continuous =
        price_floating_lookback_continuous(100.0, 90.0, 0.03, 0.0, 0.2, 1.0, OptionType::Call);
    let mut previous = f64::NEG_INFINITY;
    for observations in [4, 12, 52, 252, 100_000] {
        let value = price_floating_lookback_discrete(
            100.0,
            90.0,
            0.03,
            0.0,
            0.2,
            1.0,
            OptionType::Call,
            observations,
        );
        assert!(value + 1e-10 >= previous);
        previous = value;
    }
    assert!((continuous - previous).abs() < 2e-2);
}

#[test]
fn lookback_dual5_greeks_match_fd_including_zero_carry() {
    for &(rate, dividend) in &[(0.03, 0.01), (0.03, 0.03)] {
        let base = [100.0, 0.2, rate, dividend, 1.0];
        let floating = price_floating_lookback_continuous(
            Dual::<5>::variable(base[0], 0),
            Dual::from_f64(90.0),
            Dual::variable(base[2], 2),
            Dual::variable(base[3], 3),
            Dual::variable(base[1], 1),
            Dual::variable(base[4], 4),
            OptionType::Call,
        );
        let fixed = price_fixed_lookback_continuous(
            Dual::<5>::variable(base[0], 0),
            Dual::from_f64(110.0),
            Dual::from_f64(95.0),
            Dual::variable(base[2], 2),
            Dual::variable(base[3], 3),
            Dual::variable(base[1], 1),
            Dual::variable(base[4], 4),
            OptionType::Call,
        );
        for j in 0..5 {
            let h = 1e-5;
            let mut up = base;
            let mut down = base;
            up[j] += h;
            down[j] -= h;
            let fp = price_floating_lookback_continuous(
                up[0],
                90.0,
                up[2],
                up[3],
                up[1],
                up[4],
                OptionType::Call,
            );
            let fm = price_floating_lookback_continuous(
                down[0],
                90.0,
                down[2],
                down[3],
                down[1],
                down[4],
                OptionType::Call,
            );
            let fd_float = (fp - fm) / (2.0 * h);
            assert!((floating.d[j] - fd_float).abs() < 1e-5, "float {rate} {j}");
            let fp = price_fixed_lookback_continuous(
                up[0],
                110.0,
                95.0,
                up[2],
                up[3],
                up[1],
                up[4],
                OptionType::Call,
            );
            let fm = price_fixed_lookback_continuous(
                down[0],
                110.0,
                95.0,
                down[2],
                down[3],
                down[1],
                down[4],
                OptionType::Call,
            );
            let fd_fixed = (fp - fm) / (2.0 * h);
            assert!((fixed.d[j] - fd_fixed).abs() < 1e-5, "fixed {rate} {j}");
        }
    }
}
