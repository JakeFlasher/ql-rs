//! Eydeland & Wolyniec (2003) & Carmona–Durrleman (2003) Three-Asset 3:2:1 Refinery
//! Crack-Spread Kirk–Margrabe Moment-Matching Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct EydelandThreeAssetCrackSpreadSpec<T: Real = f64> {
    pub crude_forward_barrel_fc: T,
    pub gasoline_forward_barrel_fg: T,
    pub distillate_forward_barrel_fh: T,
    pub fixed_refining_strike_k: T,
    pub crude_vol_sigma_c: T,
    pub gasoline_vol_sigma_g: T,
    pub distillate_vol_sigma_h: T,
    pub corr_gas_dist_rho_gh: T,
    pub corr_product_crude_rho_pc: T,
    pub risk_free_rate_r: T,
    pub expiry_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct EydelandThreeAssetCrackSpreadReport<T: Real = f64> {
    pub composite_product_basket_forward: T,
    pub intrinsic_321_crack_margin: T,
    pub composite_product_volatility: T,
    pub kirk_composite_crack_volatility: T,
    pub crack_spread_call_price_per_bbl: T,
    pub crude_leg_delta_sensitivity: T,
}

pub fn evaluate_eydeland_crack_spread_three_asset_kirk<T: Real>(
    spec: &EydelandThreeAssetCrackSpreadSpec<T>,
) -> EydelandThreeAssetCrackSpreadReport<T> {
    let fc = spec.crude_forward_barrel_fc.max_val(T::from_f64(1e-4));
    let fg = spec.gasoline_forward_barrel_fg.max_val(T::from_f64(1e-4));
    let fh = spec.distillate_forward_barrel_fh.max_val(T::from_f64(1e-4));
    let k = spec.fixed_refining_strike_k;
    let t = spec.expiry_years_t.max_val(T::from_f64(0.02));
    let df = (-spec.risk_free_rate_r * t).exp();

    // 3:2:1 crack per barrel of crude: ProductBasket = (2/3) * F_gas + (1/3) * F_distillate:
    let w_g = T::from_f64(2.0 / 3.0);
    let w_h = T::from_f64(1.0 / 3.0);
    let f_prod = w_g * fg + w_h * fh;
    let intrinsic = f_prod - fc - k;

    // Levy (1992) moment-matched volatility of the 2:1 product basket:
    let share_g = (w_g * fg) / f_prod;
    let share_h = (w_h * fh) / f_prod;
    let var_prod = (share_g * share_g * spec.gasoline_vol_sigma_g * spec.gasoline_vol_sigma_g
        + share_h * share_h * spec.distillate_vol_sigma_h * spec.distillate_vol_sigma_h
        + T::from_f64(2.0)
            * share_g
            * share_h
            * spec.corr_gas_dist_rho_gh
            * spec.gasoline_vol_sigma_g
            * spec.distillate_vol_sigma_h)
        .max_val(T::from_f64(1e-8));
    let sig_prod = var_prod.sqrt();

    // Kirk (1995) approximation against composite strike A = F_crude + K:
    let a = (fc + k).max_val(T::from_f64(1e-4));
    let w_c = fc / a;
    let var_kirk = (var_prod
        - T::from_f64(2.0)
            * spec.corr_product_crude_rho_pc
            * sig_prod
            * spec.crude_vol_sigma_c
            * w_c
        + spec.crude_vol_sigma_c * spec.crude_vol_sigma_c * w_c * w_c)
        .max_val(T::from_f64(1e-8));
    let sig_kirk = var_kirk.sqrt();
    let std = sig_kirk * t.sqrt();

    let d1 = ((f_prod / a).ln() + T::from_f64(0.5) * std * std) / std;
    let d2 = d1 - std;
    let call_price = df * (f_prod * norm_cdf(d1) - a * norm_cdf(d2));
    let crude_delta = -df * norm_cdf(d2);

    EydelandThreeAssetCrackSpreadReport {
        composite_product_basket_forward: f_prod,
        intrinsic_321_crack_margin: intrinsic,
        composite_product_volatility: sig_prod,
        kirk_composite_crack_volatility: sig_kirk,
        crack_spread_call_price_per_bbl: call_price,
        crude_leg_delta_sensitivity: crude_delta,
    }
}

#[cfg(test)]
mod phase41_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn eydeland_crack_spread_three_asset_kirk_invariants_and_ad_parity() {
        let spec = EydelandThreeAssetCrackSpreadSpec {
            crude_forward_barrel_fc: 76.0,
            gasoline_forward_barrel_fg: 94.0,
            distillate_forward_barrel_fh: 98.0,
            fixed_refining_strike_k: 6.5,
            crude_vol_sigma_c: 0.34,
            gasoline_vol_sigma_g: 0.38,
            distillate_vol_sigma_h: 0.36,
            corr_gas_dist_rho_gh: 0.82,
            corr_product_crude_rho_pc: 0.86,
            risk_free_rate_r: 0.04,
            expiry_years_t: 0.5,
        };
        let rep = evaluate_eydeland_crack_spread_three_asset_kirk(&spec);
        assert!(rep.intrinsic_321_crack_margin > 0.0);
        assert!(rep.crack_spread_call_price_per_bbl > rep.intrinsic_321_crack_margin * 0.9);
        assert!(rep.crude_leg_delta_sensitivity < 0.0);

        let h = 1e-6;
        let r_ad =
            evaluate_eydeland_crack_spread_three_asset_kirk(&EydelandThreeAssetCrackSpreadSpec {
                crude_forward_barrel_fc: Dual::<1>::variable(spec.crude_forward_barrel_fc, 0),
                gasoline_forward_barrel_fg: Dual::<1>::from_value(spec.gasoline_forward_barrel_fg),
                distillate_forward_barrel_fh: Dual::<1>::from_value(
                    spec.distillate_forward_barrel_fh,
                ),
                fixed_refining_strike_k: Dual::<1>::from_value(spec.fixed_refining_strike_k),
                crude_vol_sigma_c: Dual::<1>::from_value(spec.crude_vol_sigma_c),
                gasoline_vol_sigma_g: Dual::<1>::from_value(spec.gasoline_vol_sigma_g),
                distillate_vol_sigma_h: Dual::<1>::from_value(spec.distillate_vol_sigma_h),
                corr_gas_dist_rho_gh: Dual::<1>::from_value(spec.corr_gas_dist_rho_gh),
                corr_product_crude_rho_pc: Dual::<1>::from_value(spec.corr_product_crude_rho_pc),
                risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
                expiry_years_t: Dual::<1>::from_value(spec.expiry_years_t),
            });
        let mut sp = spec;
        sp.crude_forward_barrel_fc += h;
        let mut sm = spec;
        sm.crude_forward_barrel_fc -= h;
        let fd = (evaluate_eydeland_crack_spread_three_asset_kirk(&sp)
            .crack_spread_call_price_per_bbl
            - evaluate_eydeland_crack_spread_three_asset_kirk(&sm).crack_spread_call_price_per_bbl)
            / (2.0 * h);
        assert!((r_ad.crack_spread_call_price_per_bbl.d[0] - fd).abs() < 1e-5);
    }
}
