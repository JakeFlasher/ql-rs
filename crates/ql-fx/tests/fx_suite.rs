use ql_ad::Dual;
use ql_fx::*;

fn quote() -> FxSmileQuote<f64> {
    FxSmileQuote {
        atm_vol: 0.10,
        rr_25d: 0.01,
        bf_25d: 0.005,
    }
}

#[test]
fn sticky_delta_strikes_and_vanna_volga_pillar_interpolation() {
    let sm = FxVannaVolgaSmile::new(100.0, 0.03, 0.01, 1.0, quote());
    let dp = garman_kohlhagen_spot_delta(100.0, sm.k1, sm.sigma1, 1.0, 0.03, 0.01, false);
    let dc = garman_kohlhagen_spot_delta(100.0, sm.k3, sm.sigma3, 1.0, 0.03, 0.01, true);
    assert!((dp + 0.25).abs() < 1e-12);
    assert!((dc - 0.25).abs() < 1e-12);
    let atm_call = garman_kohlhagen_spot_delta(100.0, sm.k2, sm.sigma2, 1.0, 0.03, 0.01, true);
    let atm_put = garman_kohlhagen_spot_delta(100.0, sm.k2, sm.sigma2, 1.0, 0.03, 0.01, false);
    assert!((atm_call + atm_put - (-0.0)).abs() < 1e-12);
    assert!((sm.implied_vol(sm.k1) - sm.sigma1).abs() < 1e-12);
    assert!((sm.implied_vol(sm.k2) - sm.sigma2).abs() < 1e-12);
    assert!((sm.implied_vol(sm.k3) - sm.sigma3).abs() < 1e-12);
}

#[test]
fn vanna_volga_dual3_quote_sensitivities_match_fd() {
    let base = quote();
    let strikes = [85.0, 95.0, 100.0, 110.0, 125.0];
    let h = 1e-6;
    for &k in &strikes {
        let q = FxSmileQuote {
            atm_vol: Dual::<3>::variable(base.atm_vol, 0),
            rr_25d: Dual::<3>::variable(base.rr_25d, 1),
            bf_25d: Dual::<3>::variable(base.bf_25d, 2),
        };
        let sm = FxVannaVolgaSmile::new(
            Dual::from_value(100.0),
            Dual::from_value(0.03),
            Dual::from_value(0.01),
            Dual::from_value(1.0),
            q,
        );
        let iv = sm.implied_vol(Dual::from_value(k));
        let pv = sm.price_vanna_volga(Dual::from_value(k), true);
        for j in 0..3 {
            let mut up = base;
            let mut dn = base;
            match j {
                0 => {
                    up.atm_vol += h;
                    dn.atm_vol -= h;
                }
                1 => {
                    up.rr_25d += h;
                    dn.rr_25d -= h;
                }
                _ => {
                    up.bf_25d += h;
                    dn.bf_25d -= h;
                }
            }
            let sup = FxVannaVolgaSmile::new(100.0, 0.03, 0.01, 1.0, up);
            let sdn = FxVannaVolgaSmile::new(100.0, 0.03, 0.01, 1.0, dn);
            let fd_iv = (sup.implied_vol(k) - sdn.implied_vol(k)) / (2.0 * h);
            let fd_pv =
                (sup.price_vanna_volga(k, true) - sdn.price_vanna_volga(k, true)) / (2.0 * h);
            assert!((iv.d[j] - fd_iv).abs() < 1e-6, "iv k={k} j={j}");
            assert!((pv.d[j] - fd_pv).abs() < 1e-6, "pv k={k} j={j}");
        }
    }
}
