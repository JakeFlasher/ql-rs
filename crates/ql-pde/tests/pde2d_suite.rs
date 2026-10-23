use ql_ad::Dual;
use ql_instruments::{Fd1dBlackScholesEngine, OptionType};
use ql_models::*;
use ql_pde::*;

#[test]
fn fd2d_heston_mcs_and_hv_match_cos_and_1d_american() {
    let p = HestonParams {
        v0: 0.04,
        kappa: 1.5,
        theta: 0.04,
        sigma_v: 0.3,
        rho: -0.7,
    };
    let cos = CosHestonEngine::default();
    for scheme in [AdiScheme::ModifiedCraigSneyd, AdiScheme::HundsdorferVerwer] {
        let pde = Fd2dHestonEngine {
            nx: 70,
            nv: 35,
            nt: 50,
            theta: 0.5,
            scheme,
        };
        for &strike in &[90.0, 100.0, 110.0] {
            let p_pde = pde.price_european(&p, 100.0, strike, 0.03, 0.0, 0.5, true);
            let p_cos = cos.price(&p, 100.0, strike, 0.03, 0.0, 0.5, true);
            assert!(
                (p_pde - p_cos).abs() < 3.5e-2,
                "scheme={scheme:?} strike={strike}: pde={p_pde} cos={p_cos}"
            );
        }
    }

    let pde = Fd2dHestonEngine {
        nx: 70,
        nv: 30,
        nt: 50,
        theta: 0.5,
        scheme: AdiScheme::ModifiedCraigSneyd,
    };
    let eur_put = pde.price_european(&p, 100.0, 105.0, 0.05, 0.0, 1.0, false);
    let am_put = pde.price_american(&p, 100.0, 105.0, 0.05, 0.0, 1.0, false);
    assert!(am_put > eur_put + 0.1, "am_put={am_put} eur_put={eur_put}");

    let p_bs_limit = HestonParams {
        v0: 0.04,
        kappa: 2.0,
        theta: 0.04,
        sigma_v: 1e-3,
        rho: 0.0,
    };
    let am_2d = pde.price_american(&p_bs_limit, 100.0, 100.0, 0.05, 0.0, 1.0, false);
    let am_1d = Fd1dBlackScholesEngine::new(200, 200).price_american_ad(
        100.0,
        100.0,
        1.0,
        0.20,
        0.05,
        OptionType::Put,
    );
    assert!(
        (am_2d - am_1d).abs() < 3.5e-2,
        "am_2d={am_2d} am_1d={am_1d}"
    );

    let p_dual = HestonParams {
        v0: Dual::<1>::from_value(0.04),
        kappa: Dual::<1>::from_value(1.5),
        theta: Dual::<1>::from_value(0.04),
        sigma_v: Dual::<1>::from_value(0.3),
        rho: Dual::<1>::from_value(-0.7),
    };
    let spot_ad = Dual::<1>::variable(100.0, 0);
    let val_ad = pde.price_european(
        &p_dual,
        spot_ad,
        Dual::from_value(100.0),
        Dual::from_value(0.03),
        Dual::from_value(0.0),
        Dual::from_value(0.5),
        true,
    );
    let h = 1e-4;
    let vp = pde.price_european(&p, 100.0 + h, 100.0, 0.03, 0.0, 0.5, true);
    let vm = pde.price_european(&p, 100.0 - h, 100.0, 0.03, 0.0, 0.5, true);
    let delta_fd = (vp - vm) / (2.0 * h);
    assert!(
        (val_ad.d[0] - delta_fd).abs() < 2e-3,
        "ad={} fd={delta_fd}",
        val_ad.d[0]
    );
}
