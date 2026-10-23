use ql_ad::*;
use ql_scripting::*;
use ql_vol::black_scholes_price;
#[test]
fn asian_call_single_step_matches_black_scholes() {
    let p = ScriptedProduct::asian_call(vec![1.0], 100.0, 1.0);
    let b = BytecodeProgram::compile(&p);
    assert_eq!(b.times.len(), 1);
    let e = ScriptMcEngine {
        n_paths: 32_768,
        seed: 42,
    };
    let x = e.price_gbm::<f64>(&b, 100.0, 0.2, 0.03, 0.0);
    let reference = black_scholes_price(100.0, 100.0, 0.2, 1.0, 0.03, 0.0, true);
    assert!((x - reference).abs() < 2e-2, "mc={x}, bs={reference}");
}
#[test]
fn autocallable_fuzzy_aad_matches_fd() {
    let p = BytecodeProgram::compile(&ScriptedProduct::autocallable(
        vec![0.5, 1.0],
        105.0,
        95.0,
        0.04,
        100.0,
        85.0,
        2.0,
    ));
    assert_eq!(p.times.len(), 2);
    let engine = ScriptMcEngine {
        n_paths: 2_048,
        seed: 7,
    };
    let price = engine.price_gbm(
        &p,
        Dual::<2>::variable(100.0, 0),
        Dual::<2>::variable(0.2, 1),
        Dual::from_value(0.03),
        Dual::from_value(0.0),
    );
    let h = 1e-4;
    let spot_up = engine.price_gbm(&p, 100.0 + h, 0.2, 0.03, 0.0);
    let spot_dn = engine.price_gbm(&p, 100.0 - h, 0.2, 0.03, 0.0);
    let vol_up = engine.price_gbm(&p, 100.0, 0.2 + h, 0.03, 0.0);
    let vol_dn = engine.price_gbm(&p, 100.0, 0.2 - h, 0.03, 0.0);
    assert!((price.d[0] - (spot_up - spot_dn) / (2.0 * h)).abs() < 1e-4);
    assert!((price.d[1] - (vol_up - vol_dn) / (2.0 * h)).abs() < 1e-4);
}
#[test]
fn tarn_monotone() {
    let e = ScriptMcEngine {
        n_paths: 256,
        seed: 4,
    };
    let a = e.price_gbm::<f64>(
        &BytecodeProgram::compile(&ScriptedProduct::tarn(vec![0.5, 1.], 100., 5., 1., 1.)),
        100.,
        0.2,
        0.03,
        0.,
    );
    let b = e.price_gbm::<f64>(
        &BytecodeProgram::compile(&ScriptedProduct::tarn(vec![0.5, 1.], 100., 15., 1., 1.)),
        100.,
        0.2,
        0.03,
        0.,
    );
    assert!(b >= a);
}

#[test]
fn tarn_target_knocks_out_at_the_cap() {
    let p = BytecodeProgram::compile(&ScriptedProduct::tarn(
        vec![0.5, 1.0, 1.5],
        100.0,
        5.0,
        1.0,
        0.5,
    ));
    let mut vars = vec![0.0; p.num_vars];
    let mut stack = Vec::new();
    let pv = p.eval_path(&[110.0, 110.0, 110.0], 0.0, &mut vars, &mut stack);
    assert!((pv - 6.0).abs() < 1e-12, "pv={pv}");
}
