use crate::{CosHestonEngine, HestonParams};
use nalgebra::{DMatrix, DVector};
use ql_ad::{Dual, Dual2, IftSystemSensitivity};
use ql_math::{MathError, levenberg_marquardt};
#[derive(Copy, Clone, Debug)]
pub struct HestonMarketQuote {
    pub strike: f64,
    pub expiry: f64,
    pub market_price: f64,
    pub is_call: bool,
}
#[derive(Clone, Debug)]
pub struct CalibratedHestonModel {
    pub params: HestonParams<f64>,
    pub rmse: f64,
    pub sensitivity: IftSystemSensitivity,
}
pub fn calibrate_heston(
    spot: f64,
    rate: f64,
    div: f64,
    initial: HestonParams<f64>,
    quotes: &[HestonMarketQuote],
) -> Result<CalibratedHestonModel, MathError> {
    let engine = CosHestonEngine::default();
    let x0 = DVector::from_vec(vec![
        initial.v0,
        initial.kappa,
        initial.theta,
        initial.sigma_v,
        initial.rho,
    ]);
    let f = |x: &DVector<f64>| {
        let mut d = vec![0.; quotes.len()];
        let mut j = DMatrix::zeros(quotes.len(), 5);
        for (i, q) in quotes.iter().enumerate() {
            let vars = [
                Dual::<5>::variable(x[0].max(1e-5), 0),
                Dual::variable(x[1].max(1e-5), 1),
                Dual::variable(x[2].max(1e-5), 2),
                Dual::variable(x[3].max(1e-5), 3),
                Dual::variable(x[4].clamp(-0.995, 0.995), 4),
            ];
            let p = engine.price(
                &HestonParams {
                    v0: vars[0],
                    kappa: vars[1],
                    theta: vars[2],
                    sigma_v: vars[3],
                    rho: vars[4],
                },
                Dual::from_value(spot),
                Dual::from_value(q.strike),
                Dual::from_value(rate),
                Dual::from_value(div),
                Dual::from_value(q.expiry),
                q.is_call,
            );
            d[i] = p.val - q.market_price;
            for k in 0..5 {
                j[(i, k)] = p.d[k];
            }
        }
        (DVector::from_vec(d), j)
    };
    let sol = levenberg_marquardt(x0, f, 1e-10, 100)?;
    let params = HestonParams {
        v0: sol[0].max(1e-5),
        kappa: sol[1].max(1e-5),
        theta: sol[2].max(1e-5),
        sigma_v: sol[3].max(1e-5),
        rho: sol[4].clamp(-0.995, 0.995),
    };
    let (r, j) = f(&sol);
    let rmse = (r.norm_squared() / quotes.len().max(1) as f64).sqrt();
    let mut h_opt = j.transpose() * &j;
    for q in quotes {
        let vars = [
            Dual2::<5>::variable(params.v0, 0),
            Dual2::<5>::variable(params.kappa, 1),
            Dual2::<5>::variable(params.theta, 2),
            Dual2::<5>::variable(params.sigma_v, 3),
            Dual2::<5>::variable(params.rho, 4),
        ];
        let p2 = engine.price(
            &HestonParams {
                v0: vars[0],
                kappa: vars[1],
                theta: vars[2],
                sigma_v: vars[3],
                rho: vars[4],
            },
            Dual2::from_value(spot),
            Dual2::from_value(q.strike),
            Dual2::from_value(rate),
            Dual2::from_value(div),
            Dual2::from_value(q.expiry),
            q.is_call,
        );
        let residual = p2.val - q.market_price;
        for i in 0..5 {
            for jcol in 0..5 {
                h_opt[(i, jcol)] += residual * p2.d2[i][jcol];
            }
        }
    }
    if h_opt.determinant().abs() < 1e-14 {
        for i in 0..5 {
            h_opt[(i, i)] += 1e-8;
        }
    }
    let sens = IftSystemSensitivity::new(h_opt, -j.transpose());
    Ok(CalibratedHestonModel {
        params,
        rmse,
        sensitivity: sens,
    })
}
