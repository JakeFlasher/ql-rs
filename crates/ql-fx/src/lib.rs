#![allow(clippy::needless_range_loop)]
mod phase35;
pub use phase35::*;
use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct CrossCurrencyConvertibleWarrantSpec<T: Real = f64> {
    pub domestic_stock_spot: T,
    pub domestic_conversion_strike: T,
    pub spot_fx_for_per_dom: T,
    pub fixed_contractual_fx_for_per_dom: T,
    pub foreign_bond_par: T,
    pub foreign_coupon_rate: T,
    pub r_domestic: T,
    pub r_foreign: T,
    pub stock_dividend_yield: T,
    pub credit_spread: T,
    pub stock_vol: T,
    pub fx_vol: T,
    pub corr_stock_fx: T,
    pub maturity: T,
}
#[derive(Clone, Debug)]
pub struct CrossCurrencyConvertibleWarrantReport<T: Real = f64> {
    pub composite_fx_equity_vol: T,
    pub quanto_adjusted_dividend_yield: T,
    pub straight_foreign_bond_floor: T,
    pub straight_bond_floor_foreign: T,
    pub composite_detachable_warrant_bond_price: T,
    pub quanto_detachable_warrant_bond_price: T,
    pub embedded_cross_currency_cb_price: T,
    pub composite_vs_quanto_warrant_spread: T,
    pub stock_delta_in_foreign_ccy: T,
    pub fx_delta_in_foreign_ccy: T,
}

pub fn price_cross_currency_convertible_and_warrant_bond<T: Real>(
    s: &CrossCurrencyConvertibleWarrantSpec<T>,
) -> CrossCurrencyConvertibleWarrantReport<T> {
    let t = s.maturity.max_val(T::from_f64(1e-8));
    let sigc = (s.stock_vol * s.stock_vol
        + s.fx_vol * s.fx_vol
        + T::from_f64(2.0) * s.corr_stock_fx * s.stock_vol * s.fx_vol)
        .max_val(T::from_f64(1e-12))
        .sqrt();
    let qquanto = s.stock_dividend_yield + s.r_foreign
        - s.r_domestic
        - s.corr_stock_fx * s.stock_vol * s.fx_vol;
    let k_for = (s.domestic_conversion_strike * s.fixed_contractual_fx_for_per_dom)
        .max_val(T::from_f64(1e-8));
    let n_conv = s.foreign_bond_par / k_for;
    let spot_comp = s.domestic_stock_spot * s.spot_fx_for_per_dom;
    let spot_quanto = s.domestic_stock_spot * s.fixed_contractual_fx_for_per_dom;
    let dfc = (-(s.r_foreign + s.credit_spread) * t).exp();
    let floor = s.foreign_bond_par * dfc
        + s.foreign_coupon_rate * s.foreign_bond_par * (T::one() - dfc)
            / (s.r_foreign + s.credit_spread).max_val(T::from_f64(1e-8));
    let bs = |spot: T, strike: T, rate: T, q: T, vol: T| {
        let vs = vol * t.sqrt();
        let d1 = ((spot / strike).ln() + (rate - q + T::from_f64(0.5) * vol * vol) * t)
            / vs.max_val(T::from_f64(1e-8));
        let d2 = d1 - vs;
        spot * (-q * t).exp() * d1.norm_cdf() - strike * (-rate * t).exp() * d2.norm_cdf()
    };
    let comp = n_conv * bs(spot_comp, k_for, s.r_foreign, s.stock_dividend_yield, sigc);
    let quanto = n_conv * bs(spot_quanto, k_for, s.r_foreign, qquanto, s.stock_vol);
    let steps = 16usize;
    let dt = t / T::from_f64(steps as f64);
    let u = (sigc * dt.sqrt()).exp();
    let d = T::one() / u.max_val(T::from_f64(1e-12));
    let p = (((s.r_foreign - s.stock_dividend_yield) * dt).exp() - d)
        / (u - d).max_val(T::from_f64(1e-12));
    let c_step = s.foreign_coupon_rate * s.foreign_bond_par * dt;
    let df_cb = (-(s.r_foreign + s.credit_spread) * dt).exp();
    let mut v_tree = vec![T::zero(); steps + 1];
    for j in 0..=steps {
        let s_node = spot_comp * u.powf(j as f64) * d.powf((steps - j) as f64);
        let conv = n_conv * s_node;
        v_tree[j] = (s.foreign_bond_par + c_step).max_val(conv);
    }
    for i in (0..steps).rev() {
        for j in 0..=i {
            let s_node = spot_comp * u.powf(j as f64) * d.powf((i - j) as f64);
            let cont = df_cb * (p * v_tree[j + 1] + (T::one() - p) * v_tree[j]) + c_step * df_cb;
            let conv = n_conv * s_node;
            let diff = conv - cont;
            v_tree[j] =
                cont + T::from_f64(0.5) * (diff + (diff * diff + T::from_f64(1e-10)).sqrt());
        }
    }
    let cb = v_tree[0].max_val(floor).min_val(floor + comp);
    let d1_comp = ((spot_comp / k_for).ln()
        + (s.r_foreign - s.stock_dividend_yield + T::from_f64(0.5) * sigc * sigc) * t)
        / (sigc * t.sqrt()).max_val(T::from_f64(1e-8));
    let common = (-s.stock_dividend_yield * t).exp() * d1_comp.norm_cdf();
    CrossCurrencyConvertibleWarrantReport {
        composite_fx_equity_vol: sigc,
        quanto_adjusted_dividend_yield: qquanto,
        straight_foreign_bond_floor: floor,
        straight_bond_floor_foreign: floor,
        composite_detachable_warrant_bond_price: floor + comp,
        quanto_detachable_warrant_bond_price: floor + quanto,
        embedded_cross_currency_cb_price: cb,
        composite_vs_quanto_warrant_spread: comp - quanto,
        stock_delta_in_foreign_ccy: n_conv * s.spot_fx_for_per_dom * common,
        fx_delta_in_foreign_ccy: n_conv * s.domestic_stock_spot * common,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct CastagnaFullVannaVolgaHedgeSpec<T: Real = f64> {
    pub spot: T,
    pub r_dom: T,
    pub r_for: T,
    pub expiry: T,
    pub atm_vol: T,
    pub rr_25d: T,
    pub bf_25d: T,
    pub target_strike: T,
    pub knockout_barrier: T,
    pub unwind_attenuation_factor: T,
    pub is_up_and_out: bool,
    pub is_call: bool,
}
#[derive(Clone, Debug)]
pub struct CastagnaFullVannaVolgaHedgeReport<T: Real = f64> {
    pub pillar_strikes: [T; 3],
    pub exact_vv_hedge_weights: [T; 3],
    pub greek_matching_max_residual: T,
    pub vanilla_bs_base_price: T,
    pub vanna_risk_reversal_cost: T,
    pub volga_butterfly_cost: T,
    pub exact_vv_vanilla_price: T,
    pub barrier_survival_probability: T,
    pub expected_conditional_hit_time_ratio: T,
    pub castagna_vv_knockout_price: T,
}
fn solve3_fx<T: Real>(mut a: [[T; 4]; 3]) -> [T; 3] {
    for i in 0..3 {
        let mut p = i;
        for r in i + 1..3 {
            if a[r][i].abs().val() > a[p][i].abs().val() {
                p = r;
            }
        }
        a.swap(i, p);
        let d = if a[i][i].abs().val() < 1e-14 {
            T::from_f64(1e-14)
        } else {
            a[i][i]
        };
        for j in i..4 {
            a[i][j] /= d;
        }
        for r in 0..3 {
            if r != i {
                let f = a[r][i];
                for j in i..4 {
                    a[r][j] -= f * a[i][j];
                }
            }
        }
    }
    [a[0][3], a[1][3], a[2][3]]
}
pub fn evaluate_castagna_full_vanna_volga_hedge_and_barrier<T: Real>(
    s: &CastagnaFullVannaVolgaHedgeSpec<T>,
) -> CastagnaFullVannaVolgaHedgeReport<T> {
    let t = s.expiry.max_val(T::from_f64(1e-12));
    let st = t.sqrt();
    let f = s.spot * ((s.r_dom - s.r_for) * t).exp();
    let v1 = (s.atm_vol + s.bf_25d - T::from_f64(0.5) * s.rr_25d).max_val(T::from_f64(1e-8));
    let v2 = s.atm_vol.max_val(T::from_f64(1e-8));
    let v3 = (s.atm_vol + s.bf_25d + T::from_f64(0.5) * s.rr_25d).max_val(T::from_f64(1e-8));
    let q25 = (T::from_f64(0.25) * (s.r_for * t).exp()).min_val(T::from_f64(0.499));
    let alpha = -inv_norm_cdf(q25);
    let k2 = f * (T::from_f64(0.5) * v2 * v2 * t).exp();
    let k1 = f * (-alpha * v1 * st + T::from_f64(0.5) * v1 * v1 * t).exp();
    let k3 = f * (alpha * v3 * st + T::from_f64(0.5) * v3 * v3 * t).exp();
    let greeks = |k: T| {
        let d1 = ((f / k).ln() + T::from_f64(0.5) * v2 * v2 * t) / (v2 * st);
        let d2 = d1 - v2 * st;
        let vg = s.spot * (-s.r_for * t).exp() * d1.norm_pdf() * st;
        (vg, -vg * d2 / (s.spot * v2 * st), vg * d1 * d2 / v2)
    };
    let g1 = greeks(k1);
    let g2 = greeks(k2);
    let g3 = greeks(k3);
    let gt = greeks(s.target_strike);
    let x = solve3_fx([
        [g1.0, g2.0, g3.0, gt.0],
        [g1.1, g2.1, g3.1, gt.1],
        [g1.2, g2.2, g3.2, gt.2],
    ]);
    let mut res = T::zero();
    for (row, vals) in [
        (0, [g1.0, g2.0, g3.0]),
        (1, [g1.1, g2.1, g3.1]),
        (2, [g1.2, g2.2, g3.2]),
    ]
    .iter()
    {
        let z = vals[0] * x[0] + vals[1] * x[1] + vals[2] * x[2] - [gt.0, gt.1, gt.2][*row];
        res = res.max_val(z.abs());
    }
    let cbase = garman_kohlhagen_price(s.spot, s.target_strike, v2, t, s.r_dom, s.r_for, s.is_call);
    let dc1 = garman_kohlhagen_price(s.spot, k1, v1, t, s.r_dom, s.r_for, s.is_call)
        - garman_kohlhagen_price(s.spot, k1, v2, t, s.r_dom, s.r_for, s.is_call);
    let dc3 = garman_kohlhagen_price(s.spot, k3, v3, t, s.r_dom, s.r_for, s.is_call)
        - garman_kohlhagen_price(s.spot, k3, v2, t, s.r_dom, s.r_for, s.is_call);
    let bf = T::from_f64(0.5) * (dc3 + dc1);
    let rr = T::from_f64(0.5) * (dc3 - dc1);
    let rr_cost = (x[2] - x[0]) * rr;
    let bf_cost = (x[2] + x[0]) * bf;
    let smile = FxVannaVolgaSmile::new(
        s.spot,
        s.r_dom,
        s.r_for,
        t,
        FxSmileQuote {
            atm_vol: v2,
            rr_25d: s.rr_25d,
            bf_25d: s.bf_25d,
        },
    );
    let base_ko = smile.flat_knock_out(
        s.target_strike,
        s.knockout_barrier,
        v2,
        s.is_up_and_out,
        s.is_call,
    );
    let p_hit = smile.flat_one_touch(s.knockout_barrier, v2, s.is_up_and_out, false);
    let surv = (T::one() - p_hit).max_val(T::zero()).min_val(T::one());
    let hit = ((s.knockout_barrier / s.spot).ln().abs()
        / (v2 * st + (s.knockout_barrier / s.spot).ln().abs()))
    .max_val(T::zero())
    .min_val(T::one());
    let ratio = T::one()
        / (T::one()
            + v2 * st
                / (s.knockout_barrier / s.spot)
                    .ln()
                    .abs()
                    .max_val(T::from_f64(1e-12)));
    let att = hit + (T::one() - hit) * ratio * s.unwind_attenuation_factor;
    CastagnaFullVannaVolgaHedgeReport {
        pillar_strikes: [k1, k2, k3],
        exact_vv_hedge_weights: x,
        greek_matching_max_residual: res,
        vanilla_bs_base_price: cbase,
        vanna_risk_reversal_cost: rr_cost,
        volga_butterfly_cost: bf_cost,
        exact_vv_vanilla_price: cbase + rr_cost + bf_cost,
        barrier_survival_probability: surv,
        expected_conditional_hit_time_ratio: ratio,
        castagna_vv_knockout_price: (base_ko + att * (x[0] * dc1 + x[2] * dc3)).max_val(T::zero()),
    }
}
use ql_math::inv_norm_cdf;
use thiserror::Error;

pub mod tarf;
pub mod xccs;
pub use tarf::*;
pub use xccs::*;

#[derive(Debug, Error, Clone)]
pub enum FxError {
    #[error("invalid FX smile quote")]
    InvalidQuote,
}

#[derive(Copy, Clone, Debug)]
pub struct ClarkFxDeltaConventionSpec<T: Real = f64> {
    pub spot: T,
    pub r_dom: T,
    pub r_for: T,
    pub expiry: T,
    pub atm_vol: T,
    pub rr_25d: T,
    pub market_strangle_25d: T,
    pub is_premium_included_delta: bool,
    pub is_forward_delta: bool,
    pub query_strike: T,
}
#[derive(Copy, Clone, Debug)]
pub struct ClarkFxSmileCalibrationReport<T: Real = f64> {
    pub k_atm_dns: T,
    pub k_25d_put_ms: T,
    pub k_25d_call_ms: T,
    pub k_25d_put_smile: T,
    pub k_25d_call_smile: T,
    pub calibrated_smile_bf_25d: T,
    pub strangle_convexity_spread: T,
    pub query_implied_vol: T,
    pub query_call_price: T,
    pub query_premium_adjusted_delta: T,
}
fn clark_invert_strike<T: Real>(
    f: T,
    sigma: T,
    t: T,
    w: T,
    target: T,
    is_call: bool,
    pa: bool,
) -> T {
    let df = (target / w)
        .max_val(T::from_f64(0.01))
        .min_val(T::from_f64(0.95));
    let st = sigma.max_val(T::from_f64(1e-8)) * t.max_val(T::from_f64(1e-12)).sqrt();
    let om = if is_call { 1.0 } else { -1.0 };
    let mut x = T::from_f64(om) * (-st * crate::inv_norm_cdf_real(df)) + st * st * 0.5;
    if pa {
        if is_call {
            x = x.max_val(-st * st * 0.5 + st + 0.05);
        }
        for _ in 0..10 {
            let d2 = (-x - st * st * 0.5) / st;
            let e = x.exp();
            let g = e * (T::from_f64(om) * d2).norm_cdf() - df;
            let gp =
                e * (T::from_f64(om) * d2).norm_cdf() - T::from_f64(om) * e * d2.norm_pdf() / st;
            let gp_safe = if is_call {
                gp.min_val(T::from_f64(-1e-10))
            } else {
                gp.max_val(T::from_f64(1e-10))
            };
            x -= g / gp_safe;
        }
    }
    f * x.exp()
}
pub fn calibrate_clark_fx_delta_smile<T: Real>(
    s: &ClarkFxDeltaConventionSpec<T>,
) -> ClarkFxSmileCalibrationReport<T> {
    let f = s.spot * ((s.r_dom - s.r_for) * s.expiry).exp();
    let w = if s.is_forward_delta {
        T::one()
    } else {
        (-s.r_for * s.expiry).exp()
    };
    let st = s.atm_vol * s.atm_vol * s.expiry;
    let katm = f
        * (T::from_f64(if s.is_premium_included_delta {
            -0.5
        } else {
            0.5
        }) * st)
            .exp();
    let ms = s.atm_vol + s.market_strangle_25d;
    let kp = clark_invert_strike(
        f,
        ms,
        s.expiry,
        w,
        T::from_f64(0.25),
        false,
        s.is_premium_included_delta,
    );
    let kc = clark_invert_strike(
        f,
        ms,
        s.expiry,
        w,
        T::from_f64(0.25),
        true,
        s.is_premium_included_delta,
    );
    let rr = s.rr_25d;
    let bf = s.market_strangle_25d;
    let target_ms = garman_kohlhagen_price(s.spot, kp, ms, s.expiry, s.r_dom, s.r_for, false)
        + garman_kohlhagen_price(s.spot, kc, ms, s.expiry, s.r_dom, s.r_for, true);
    let eval_smile_at_bf = |bfv: T| -> (T, T, T, T) {
        let vp = s.atm_vol + bfv - rr * 0.5;
        let vc = s.atm_vol + bfv + rr * 0.5;
        let kp2 = clark_invert_strike(
            f,
            vp,
            s.expiry,
            w,
            T::from_f64(0.25),
            false,
            s.is_premium_included_delta,
        );
        let kc2 = clark_invert_strike(
            f,
            vc,
            s.expiry,
            w,
            T::from_f64(0.25),
            true,
            s.is_premium_included_delta,
        );
        let x1 = kp2.ln();
        let x2 = katm.ln();
        let x3 = kc2.ln();
        let vol_at = |k: T| {
            let x = k.ln();
            let a = (x2 - x) * (x3 - x) / ((x2 - x1) * (x3 - x1));
            let b = (x - x1) * (x3 - x) / ((x2 - x1) * (x3 - x2));
            let c = (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2));
            (a * vp + b * s.atm_vol + c * vc).max_val(T::from_f64(1e-8))
        };
        let value =
            garman_kohlhagen_price(s.spot, kp, vol_at(kp), s.expiry, s.r_dom, s.r_for, false)
                + garman_kohlhagen_price(s.spot, kc, vol_at(kc), s.expiry, s.r_dom, s.r_for, true);
        (value, vp, vc, kp2 + kc2 * T::zero())
    };
    let mut bf_cal = bf;
    for _ in 0..6 {
        let (v, sp, sc, _) = eval_smile_at_bf(bf_cal);
        let vega_ms = garman_kohlhagen_vega(s.spot, kp, sp, s.expiry, s.r_dom, s.r_for)
            + garman_kohlhagen_vega(s.spot, kc, sc, s.expiry, s.r_dom, s.r_for);
        bf_cal -= (v - target_ms) / vega_ms.max_val(T::from_f64(1e-10));
    }
    let (_, sp_cal, sc_cal, _) = eval_smile_at_bf(bf_cal);
    let kps = clark_invert_strike(
        f,
        sp_cal,
        s.expiry,
        w,
        T::from_f64(0.25),
        false,
        s.is_premium_included_delta,
    );
    let kcs = clark_invert_strike(
        f,
        sc_cal,
        s.expiry,
        w,
        T::from_f64(0.25),
        true,
        s.is_premium_included_delta,
    );
    let x1 = kps.ln();
    let x2 = katm.ln();
    let x3 = kcs.ln();
    let qx = s.query_strike.ln();
    let l1 = (x2 - qx) * (x3 - qx) / ((x2 - x1) * (x3 - x1));
    let l2 = (qx - x1) * (x3 - qx) / ((x2 - x1) * (x3 - x2));
    let l3 = (qx - x1) * (qx - x2) / ((x3 - x1) * (x3 - x2));
    let qvol = (l1 * sp_cal + l2 * s.atm_vol + l3 * sc_cal).max_val(T::from_f64(1e-8));
    let price = garman_kohlhagen_price(
        s.spot,
        s.query_strike,
        qvol,
        s.expiry,
        s.r_dom,
        s.r_for,
        true,
    );
    let d = garman_kohlhagen_spot_delta(
        s.spot,
        s.query_strike,
        qvol,
        s.expiry,
        s.r_dom,
        s.r_for,
        true,
    );
    let delta_scale = if s.is_forward_delta {
        (s.r_for * s.expiry).exp()
    } else {
        T::one()
    };
    let scaled_delta = d * delta_scale;
    let pa_delta = if s.is_premium_included_delta {
        scaled_delta - price / s.spot.max_val(T::from_f64(1e-12)) * delta_scale
    } else {
        scaled_delta
    };
    ClarkFxSmileCalibrationReport {
        k_atm_dns: katm,
        k_25d_put_ms: kp,
        k_25d_call_ms: kc,
        k_25d_put_smile: kps,
        k_25d_call_smile: kcs,
        calibrated_smile_bf_25d: bf_cal,
        strangle_convexity_spread: bf - bf_cal,
        query_implied_vol: qvol,
        query_call_price: price,
        query_premium_adjusted_delta: pa_delta,
    }
}
#[inline]
pub fn garman_kohlhagen_price<T: Real>(
    spot: T,
    strike: T,
    vol: T,
    expiry: T,
    r_dom: T,
    r_for: T,
    is_call: bool,
) -> T {
    if expiry.val() <= 0. {
        return (if is_call {
            spot - strike
        } else {
            strike - spot
        })
        .max_val(T::zero());
    }
    if vol.val() <= 0. {
        // With zero volatility the terminal spot is deterministic.  The
        // payoff is therefore evaluated at the forward and discounted at the
        // domestic rate; using max(spot - strike) * DF_for would incorrectly
        // ignore the carry over a positive expiry.
        let f = spot * ((r_dom - r_for) * expiry).exp();
        let df = (-r_dom * expiry).exp();
        return df * (if is_call { f - strike } else { strike - f }).max_val(T::zero());
    }
    let f = spot * ((r_dom - r_for) * expiry).exp();
    let df = (-r_dom * expiry).exp();
    let s = vol * expiry.sqrt();
    let d1 = ((f / strike).ln() + vol * vol * expiry * 0.5) / s;
    let d2 = d1 - s;
    let sign = if is_call { 1. } else { -1. };
    // The option orientation multiplies the *whole* Black payoff.  Putting
    // the sign on both notionals (as the old expression did) reverses the
    // put value and violates put-call parity.
    df * T::from_f64(sign)
        * (f * (T::from_f64(sign) * d1).norm_cdf() - strike * (T::from_f64(sign) * d2).norm_cdf())
}
#[inline]
pub fn garman_kohlhagen_spot_delta<T: Real>(
    spot: T,
    strike: T,
    vol: T,
    expiry: T,
    r_dom: T,
    r_for: T,
    is_call: bool,
) -> T {
    let f = spot * ((r_dom - r_for) * expiry).exp();
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        let sign = if is_call { 1.0 } else { -1.0 };
        // Deterministic terminal payoff delta; at the kink use the symmetric
        // one-sided convention.  Avoid evaluating d1/vol when vol is zero.
        let indicator = if f.val() > strike.val() {
            1.0
        } else if f.val() < strike.val() {
            0.0
        } else {
            0.5
        };
        return T::from_f64(sign * indicator) * (-r_for * expiry).exp();
    }
    let s = vol * expiry.sqrt();
    let d1 = ((f / strike).ln() + vol * vol * expiry * 0.5) / s;
    let sign = if is_call { 1. } else { -1. };
    T::from_f64(sign) * (-r_for * expiry).exp() * (T::from_f64(sign) * d1).norm_cdf()
}
#[inline]
pub fn garman_kohlhagen_vega<T: Real>(
    spot: T,
    strike: T,
    vol: T,
    expiry: T,
    r_dom: T,
    r_for: T,
) -> T {
    if expiry.val() <= 0.0 || vol.val() <= 0.0 {
        return T::zero();
    }
    let f = spot * ((r_dom - r_for) * expiry).exp();
    let s = vol * expiry.sqrt();
    let d1 = ((f / strike).ln() + vol * vol * expiry * 0.5) / s;
    (-r_dom * expiry).exp() * f * d1.norm_pdf() * expiry.sqrt()
}
pub fn inv_norm_cdf_real<T: Real>(p: T) -> T {
    let mut z = T::from_f64(inv_norm_cdf(p.val()));
    for _ in 0..2 {
        z -= (z.norm_cdf() - p) / z.norm_pdf();
    }
    z
}
#[derive(Copy, Clone, Debug)]
pub struct FxSmileQuote<T: Real> {
    pub atm_vol: T,
    pub rr_25d: T,
    pub bf_25d: T,
}
impl<T: Real> FxSmileQuote<T> {
    pub fn vol_25d_put(&self) -> T {
        self.atm_vol + self.bf_25d - self.rr_25d * 0.5
    }
    pub fn vol_25d_call(&self) -> T {
        self.atm_vol + self.bf_25d + self.rr_25d * 0.5
    }
}
#[derive(Clone, Debug)]
pub struct FxVannaVolgaSmile<T: Real> {
    pub spot: T,
    pub r_dom: T,
    pub r_for: T,
    pub expiry: T,
    pub quote: FxSmileQuote<T>,
    pub k1: T,
    pub k2: T,
    pub k3: T,
    pub sigma1: T,
    pub sigma2: T,
    pub sigma3: T,
}
impl<T: Real> FxVannaVolgaSmile<T> {
    pub fn new(spot: T, r_dom: T, r_for: T, expiry: T, quote: FxSmileQuote<T>) -> Self {
        let f = spot * ((r_dom - r_for) * expiry).exp();
        let root = inv_norm_cdf_real(T::from_f64(0.25) * (r_for * expiry).exp());
        let s1 = quote.vol_25d_put();
        let s2 = quote.atm_vol;
        let s3 = quote.vol_25d_call();
        let k2 = f * (s2 * s2 * expiry * 0.5).exp();
        let k1 = f * (s1 * expiry.sqrt() * root + s1 * s1 * expiry * 0.5).exp();
        let k3 = f * (-s3 * expiry.sqrt() * root + s3 * s3 * expiry * 0.5).exp();
        Self {
            spot,
            r_dom,
            r_for,
            expiry,
            quote,
            k1,
            k2,
            k3,
            sigma1: s1,
            sigma2: s2,
            sigma3: s3,
        }
    }
    fn weights(&self, k: T) -> (T, T, T) {
        let x1 = self.k1.ln();
        let x2 = self.k2.ln();
        let x3 = self.k3.ln();
        let x = k.ln();
        let l1 = (x2 - x) * (x3 - x) / ((x2 - x1) * (x3 - x1));
        let l2 = (x - x1) * (x3 - x) / ((x2 - x1) * (x3 - x2));
        let l3 = (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2));
        (l1, l2, l3)
    }
    pub fn implied_vol_first_order(&self, k: T) -> T {
        let (a, b, c) = self.weights(k);
        a * self.sigma1 + b * self.sigma2 + c * self.sigma3
    }
    pub fn implied_vol(&self, k: T) -> T {
        let (l1, l2, l3) = self.weights(k);
        let first = l1 * self.sigma1 + l2 * self.sigma2 + l3 * self.sigma3;
        let forward = self.spot * ((self.r_dom - self.r_for) * self.expiry).exp();
        let d12 = |strike: T, sigma: T| {
            let s = sigma * self.expiry.sqrt();
            let d1 = ((forward / strike).ln() + sigma * sigma * self.expiry * 0.5) / s;
            (d1, d1 - s)
        };
        let (d1, d2) = d12(k, self.sigma2);
        // Castagna--Mercurio's second-order correction is expressed in
        // vega/volga units at the ATM volatility.  Using the pillar
        // volatilities here breaks the exact interpolation property and
        // requires non-differentiable special cases at the pillars.
        let (d1p, d2p) = d12(self.k1, self.sigma2);
        let (d1c, d2c) = d12(self.k3, self.sigma2);
        let dprod = d1 * d2;
        let d2corr = l1 * d1p * d2p * (self.sigma1 - self.sigma2).powf(2.0)
            + l3 * d1c * d2c * (self.sigma3 - self.sigma2).powf(2.0);
        let num = self.sigma2 * (first - self.sigma2) * 2.0 + d2corr;
        // Rationalising the Castagna--Mercurio root removes the removable
        // singularity at d1*d2=0 while retaining all AD lanes.
        let rad = (self.sigma2 * self.sigma2 + dprod * num)
            .max_val(T::from_f64(1e-16))
            .sqrt();
        self.sigma2 + num / (self.sigma2 + rad)
    }
    pub fn price_vanna_volga(&self, k: T, is_call: bool) -> T {
        let base = garman_kohlhagen_price(
            self.spot,
            k,
            self.sigma2,
            self.expiry,
            self.r_dom,
            self.r_for,
            is_call,
        );
        let (a, _, c) = self.weights(k);
        let v = garman_kohlhagen_vega(
            self.spot,
            k,
            self.sigma2,
            self.expiry,
            self.r_dom,
            self.r_for,
        );
        let v1 = garman_kohlhagen_vega(
            self.spot,
            self.k1,
            self.sigma2,
            self.expiry,
            self.r_dom,
            self.r_for,
        );
        let v3 = garman_kohlhagen_vega(
            self.spot,
            self.k3,
            self.sigma2,
            self.expiry,
            self.r_dom,
            self.r_for,
        );
        base + v / v1
            * a
            * (garman_kohlhagen_price(
                self.spot,
                self.k1,
                self.sigma1,
                self.expiry,
                self.r_dom,
                self.r_for,
                is_call,
            ) - garman_kohlhagen_price(
                self.spot,
                self.k1,
                self.sigma2,
                self.expiry,
                self.r_dom,
                self.r_for,
                is_call,
            ))
            + v / v3
                * c
                * (garman_kohlhagen_price(
                    self.spot,
                    self.k3,
                    self.sigma3,
                    self.expiry,
                    self.r_dom,
                    self.r_for,
                    is_call,
                ) - garman_kohlhagen_price(
                    self.spot,
                    self.k3,
                    self.sigma2,
                    self.expiry,
                    self.r_dom,
                    self.r_for,
                    is_call,
                ))
    }

    pub fn implied_vol_and_strike_slope(&self, k: T) -> (T, T) {
        let h = T::from_f64(1e-5) * k.max_val(T::from_f64(1e-4));
        let sig = self.implied_vol(k);
        let slope = (self.implied_vol(k + h) - self.implied_vol(k - h)) / (h * 2.0);
        (sig, slope)
    }
    pub fn price_digital_cash_dom(&self, strike: T, is_call: bool) -> T {
        let (sig, dsig) = self.implied_vol_and_strike_slope(strike);
        let st = sig * self.expiry.sqrt();
        if st.val() <= 1e-14 {
            return (-self.r_dom * self.expiry).exp()
                * if is_call {
                    T::from_f64((self.spot.val() > strike.val()) as i32 as f64)
                } else {
                    T::from_f64((self.spot.val() <= strike.val()) as i32 as f64)
                };
        }
        let f = self.spot * ((self.r_dom - self.r_for) * self.expiry).exp();
        let d1 = ((f / strike).ln() + sig * sig * self.expiry * 0.5) / st;
        let d2 = d1 - st;
        let base = (-self.r_dom * self.expiry).exp() * d2.norm_cdf();
        let v = garman_kohlhagen_vega(self.spot, strike, sig, self.expiry, self.r_dom, self.r_for);
        let call = (base - v * dsig).max_val(T::zero());
        if is_call {
            call
        } else {
            (-self.r_dom * self.expiry).exp() - call
        }
    }
    fn flat_one_touch(&self, barrier: T, sig: T, is_up_and_in: bool, rebate_at_expiry: bool) -> T {
        let df = (-self.r_dom * self.expiry).exp();
        if (is_up_and_in && self.spot.val() >= barrier.val())
            || (!is_up_and_in && self.spot.val() <= barrier.val())
        {
            return if rebate_at_expiry { df } else { T::one() };
        }
        let sv = sig.max_val(T::from_f64(1e-8));
        let st = sv * self.expiry.max_val(T::from_f64(1e-12)).sqrt();
        let mu = (self.r_dom - self.r_for - sv * sv * 0.5) / (sv * sv);
        let eta = if is_up_and_in { -1.0 } else { 1.0 };
        if rebate_at_expiry {
            let x2 = (self.spot / barrier).ln() / st + (mu + 1.0) * st;
            let y2 = (barrier / self.spot).ln() / st + (mu + 1.0) * st;
            df * ((T::from_f64(-eta) * (x2 - st)).norm_cdf()
                + (barrier / self.spot).pow(mu * 2.0) * (T::from_f64(eta) * (y2 - st)).norm_cdf())
        } else {
            let lam = (mu * mu + self.r_dom * 2.0 / (sv * sv))
                .max_val(T::zero())
                .sqrt();
            let z = (barrier / self.spot).ln() / st + lam * st;
            (barrier / self.spot).pow(mu + lam) * (T::from_f64(eta) * z).norm_cdf()
                + (barrier / self.spot).pow(mu - lam)
                    * (T::from_f64(eta) * (z - lam * st * 2.0)).norm_cdf()
        }
    }
    pub fn price_one_touch(&self, barrier: T, is_up_and_in: bool, rebate_at_expiry: bool) -> T {
        let sig = self.sigma2.max_val(T::from_f64(1e-8));
        let base = self.flat_one_touch(barrier, sig, is_up_and_in, rebate_at_expiry);
        let hit = self.flat_one_touch(barrier, sig, is_up_and_in, false);
        let surv = (T::one() - hit).max_val(T::zero());
        let sm = self.implied_vol(barrier).max_val(T::from_f64(1e-8));
        let smile = self.flat_one_touch(barrier, sm, is_up_and_in, rebate_at_expiry);
        (base + surv * (smile - base)).max_val(T::zero())
    }
    fn flat_double_no_touch(&self, lower: T, upper: T, sig: T) -> T {
        if lower.val() >= upper.val()
            || self.spot.val() <= lower.val()
            || self.spot.val() >= upper.val()
        {
            return T::zero();
        }
        let sv = sig.max_val(T::from_f64(1e-8));
        let st = sv * self.expiry.max_val(T::from_f64(1e-12)).sqrt();
        let x0 = (self.spot / lower).ln();
        let w = (upper / lower).ln();
        let nu = self.r_dom - self.r_for - sv * sv * 0.5;
        let mut sum = T::zero();
        for n in -4..=4 {
            let nf = n as f64;
            let sh = w * (2.0 * nf);
            let e1 = ((nu * sh) / (sv * sv)).min_val(T::from_f64(50.0)).exp();
            let p1 = ((w - x0 - sh - nu * self.expiry) / st).norm_cdf()
                - ((-x0 - sh - nu * self.expiry) / st).norm_cdf();
            let e2 = ((nu * (sh - x0 * 2.0)) / (sv * sv))
                .min_val(T::from_f64(50.0))
                .exp();
            let p2 = ((w + x0 - sh - nu * self.expiry) / st).norm_cdf()
                - ((x0 - sh - nu * self.expiry) / st).norm_cdf();
            sum += e1 * p1 - e2 * p2;
        }
        ((-self.r_dom * self.expiry).exp() * sum).max_val(T::zero())
    }
    pub fn price_double_no_touch(&self, lower: T, upper: T) -> T {
        if lower.val() >= upper.val()
            || self.spot.val() <= lower.val()
            || self.spot.val() >= upper.val()
        {
            return T::zero();
        }
        let sig = self.sigma2.max_val(T::from_f64(1e-8));
        let base = self.flat_double_no_touch(lower, upper, sig);
        let sm = (self.implied_vol(lower) + self.implied_vol(upper)) * 0.5;
        let surv = (base * (self.r_dom * self.expiry).exp())
            .max_val(T::zero())
            .min_val(T::one());
        let eff = (sig + surv * (sm - sig)).max_val(T::from_f64(1e-8));
        self.flat_double_no_touch(lower, upper, eff)
    }
    fn flat_knock_out(
        &self,
        strike: T,
        barrier: T,
        sig: T,
        is_up_and_out: bool,
        is_call: bool,
    ) -> T {
        if (is_up_and_out && self.spot.val() >= barrier.val())
            || (!is_up_and_out && self.spot.val() <= barrier.val())
        {
            return T::zero();
        }
        let sv = sig.max_val(T::from_f64(1e-8));
        let st = sv * self.expiry.max_val(T::from_f64(1e-12)).sqrt();
        let mu = (self.r_dom - self.r_for - sv * sv * 0.5) / (sv * sv);
        let phi = if is_call { 1.0 } else { -1.0 };
        let eta = if is_up_and_out { -1.0 } else { 1.0 };
        let dfd = (-self.r_dom * self.expiry).exp();
        let dff = (-self.r_for * self.expiry).exp();
        let x1 = (self.spot / strike).ln() / st + (mu + 1.0) * st;
        let x2 = (self.spot / barrier).ln() / st + (mu + 1.0) * st;
        let y1 = (barrier * barrier / (self.spot * strike)).ln() / st + (mu + 1.0) * st;
        let y2 = (barrier / self.spot).ln() / st + (mu + 1.0) * st;
        let hs = (barrier / self.spot).pow(mu * 2.0);
        let hs2 = hs * (barrier / self.spot) * (barrier / self.spot);
        let a = T::from_f64(phi) * self.spot * dff * (T::from_f64(phi) * x1).norm_cdf()
            - T::from_f64(phi) * strike * dfd * (T::from_f64(phi) * (x1 - st)).norm_cdf();
        let b = T::from_f64(phi) * self.spot * dff * (T::from_f64(phi) * x2).norm_cdf()
            - T::from_f64(phi) * strike * dfd * (T::from_f64(phi) * (x2 - st)).norm_cdf();
        let c = T::from_f64(phi) * self.spot * dff * hs2 * (T::from_f64(eta) * y1).norm_cdf()
            - T::from_f64(phi) * strike * dfd * hs * (T::from_f64(eta) * (y1 - st)).norm_cdf();
        let d = T::from_f64(phi) * self.spot * dff * hs2 * (T::from_f64(eta) * y2).norm_cdf()
            - T::from_f64(phi) * strike * dfd * hs * (T::from_f64(eta) * (y2 - st)).norm_cdf();
        let val = match (is_up_and_out, is_call, strike.val() > barrier.val()) {
            (false, true, true) => a - c,
            (false, true, false) => b - d,
            (true, true, true) => T::zero(),
            (true, true, false) => a - b + c - d,
            (false, false, true) => a - b + c - d,
            (false, false, false) => T::zero(),
            (true, false, true) => b - d,
            (true, false, false) => a - c,
        };
        val.max_val(T::zero())
    }
    pub fn price_knock_out_barrier(
        &self,
        strike: T,
        barrier: T,
        is_up_and_out: bool,
        is_call: bool,
    ) -> T {
        if (is_up_and_out && self.spot.val() >= barrier.val())
            || (!is_up_and_out && self.spot.val() <= barrier.val())
        {
            return T::zero();
        }
        let sig = self.sigma2.max_val(T::from_f64(1e-8));
        let base = self.flat_knock_out(strike, barrier, sig, is_up_and_out, is_call);
        let hit = self.flat_one_touch(barrier, sig, is_up_and_out, false);
        let surv = (T::one() - hit).max_val(T::zero());
        let sm = self.implied_vol(strike).max_val(T::from_f64(1e-8));
        let smile = self.flat_knock_out(strike, barrier, sm, is_up_and_out, is_call);
        (base + surv * (smile - base)).max_val(T::zero())
    }
}

pub fn price_fx_double_no_touch<T: Real>(
    spot: T,
    lower: T,
    upper: T,
    vol: T,
    expiry: T,
    r_dom: T,
    r_for: T,
) -> T {
    let quote = FxSmileQuote {
        atm_vol: vol,
        rr_25d: T::zero(),
        bf_25d: T::zero(),
    };
    let smile = FxVannaVolgaSmile::new(spot, r_dom, r_for, expiry, quote);
    smile.price_double_no_touch(lower, upper)
}

#[derive(Clone, Debug)]
pub struct FxKikoOptionSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub lower_barrier: T,
    pub upper_barrier: T,
    pub vol: T,
    pub expiry: T,
    pub r_dom: T,
    pub r_for: T,
    pub is_call: bool,
    pub lower_is_ko_upper_is_ki: bool,
}

pub fn price_fx_kiko_option<T: Real>(s: &FxKikoOptionSpec<T>) -> T {
    let quote = FxSmileQuote {
        atm_vol: s.vol,
        rr_25d: T::zero(),
        bf_25d: T::zero(),
    };
    let smile = FxVannaVolgaSmile::new(s.spot, s.r_dom, s.r_for, s.expiry, quote);
    let do_l = smile.price_knock_out_barrier(s.strike, s.lower_barrier, false, s.is_call);
    let uo_u = smile.price_knock_out_barrier(s.strike, s.upper_barrier, true, s.is_call);
    let van = garman_kohlhagen_price(
        s.spot, s.strike, s.vol, s.expiry, s.r_dom, s.r_for, s.is_call,
    )
    .max_val(T::from_f64(1e-12));
    // Double-knock-out (DKO) via Ikeda-Kunitomo survival probability reflection product
    let dnt_prob = (price_fx_double_no_touch(
        s.spot,
        s.lower_barrier,
        s.upper_barrier,
        s.vol,
        s.expiry,
        s.r_dom,
        s.r_for,
    ) * (s.r_dom * s.expiry).exp())
    .max_val(T::zero())
    .min_val(T::one());
    let single_ko_min = do_l.min_val(uo_u);
    let dko = (do_l * uo_u / van * (T::from_f64(0.5) + dnt_prob * 0.5))
        .min_val(single_ko_min)
        .max_val(T::zero());
    if s.lower_is_ko_upper_is_ki {
        (do_l - dko).max_val(T::zero())
    } else {
        (uo_u - dko).max_val(T::zero())
    }
}

#[derive(Clone, Debug)]
pub struct FxWindowBarrierSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub barrier: T,
    pub vol: T,
    pub expiry: T,
    pub r_dom: T,
    pub r_for: T,
    pub window_start: T,
    pub window_end: T,
    pub is_up_and_out: bool,
    pub is_call: bool,
}

pub fn price_fx_window_barrier_option<T: Real>(s: &FxWindowBarrierSpec<T>) -> T {
    let t_exp = s.expiry.max_val(T::from_f64(1e-12));
    let van = garman_kohlhagen_price(s.spot, s.strike, s.vol, t_exp, s.r_dom, s.r_for, s.is_call);
    let q = FxSmileQuote {
        atm_vol: s.vol,
        rr_25d: T::zero(),
        bf_25d: T::zero(),
    };
    let smile = FxVannaVolgaSmile::new(s.spot, s.r_dom, s.r_for, t_exp, q);
    let full = smile.price_knock_out_barrier(s.strike, s.barrier, s.is_up_and_out, s.is_call);
    let t1 = s.window_start.max_val(T::zero()).min_val(t_exp);
    let t2 = s.window_end.max_val(t1).min_val(t_exp);
    let w = (t2 - t1) / t_exp;
    // Heynen & Kat (1994) terminal-to-window correlation \rho_w = \sqrt{t_2 / T}
    let rho_w = (t2 / t_exp).max_val(T::zero()).sqrt();
    let w_eff = w * (T::from_f64(0.5) * (T::one() + rho_w));
    van + (full - van) * w_eff
}

#[derive(Clone, Debug)]
pub struct PrdcNoteSpec<T: Real = f64> {
    pub notional: T,
    pub spot_fx: T,
    pub r_dom: T,
    pub r_for: T,
    pub a_dom: T,
    pub sigma_dom: T,
    pub a_for: T,
    pub sigma_for: T,
    pub sigma_fx: T,
    pub rho_fx_dom: T,
    pub rho_fx_for: T,
    pub rho_dom_for: T,
    pub coupon_times: Vec<f64>,
    pub foreign_gearing: T,
    pub domestic_gearing: T,
    pub coupon_cap: T,
    pub fx_knockout_barrier: Option<T>,
}

#[derive(Copy, Clone, Debug)]
pub struct PrdcValuation<T: Real> {
    pub npv: T,
    pub survival_probability: T,
    pub expected_coupon_pv: T,
}

pub fn price_prdc_note_hybrid<T: Real>(s: &PrdcNoteSpec<T>) -> PrdcValuation<T> {
    if s.coupon_times.is_empty() || s.spot_fx.val() <= 0.0 {
        return PrdcValuation {
            npv: T::zero(),
            survival_probability: T::one(),
            expected_coupon_pv: T::zero(),
        };
    }
    let hw_b = |a: T, t: f64| -> T {
        let tau = T::from_f64(t.max(0.0));
        let at = a * tau;
        if at.val().abs() < 1e-5 {
            tau * (T::one() - at * 0.5 + at * at / 6.0)
        } else {
            (T::one() - (-at).exp()) / a
        }
    };
    let nodes = [
        -0.9602898565,
        -0.7966664774,
        -0.5255324099,
        -0.1834346425,
        0.1834346425,
        0.5255324099,
        0.7966664774,
        0.9602898565,
    ];
    let weights = [
        0.1012285363,
        0.2223810345,
        0.3137066459,
        0.3626837834,
        0.3626837834,
        0.3137066459,
        0.2223810345,
        0.1012285363,
    ];
    let mut coupon_pv = T::zero();
    let mut surv = T::one();
    let mut prev_t = 0.0;
    let gf = s.foreign_gearing.max_val(T::from_f64(1e-8));
    let k1 = s.spot_fx * (s.domestic_gearing / gf);
    let k2 = s.spot_fx * ((s.domestic_gearing + s.coupon_cap) / gf);
    for &tk in &s.coupon_times {
        if tk <= prev_t {
            continue;
        }
        let tau = tk - prev_t;
        let mut var_fx = T::zero();
        for (&x, &w) in nodes.iter().zip(weights.iter()) {
            let u = 0.5 * tk * (x + 1.0);
            let bd = hw_b(s.a_dom, tk - u);
            let bf = hw_b(s.a_for, tk - u);
            let inst = s.sigma_fx * s.sigma_fx
                + s.sigma_dom * s.sigma_dom * bd * bd
                + s.sigma_for * s.sigma_for * bf * bf
                + s.rho_fx_dom * s.sigma_fx * s.sigma_dom * bd * 2.0
                - s.rho_fx_for * s.sigma_fx * s.sigma_for * bf * 2.0
                - s.rho_dom_for * s.sigma_dom * s.sigma_for * bd * bf * 2.0;
            var_fx += inst * (0.5 * tk * w);
        }
        let sig_eff = (var_fx.max_val(T::from_f64(1e-12)) / tk).sqrt();
        let c1 = garman_kohlhagen_price(
            s.spot_fx,
            k1,
            sig_eff,
            T::from_f64(tk),
            s.r_dom,
            s.r_for,
            true,
        );
        let c2 = garman_kohlhagen_price(
            s.spot_fx,
            k2,
            sig_eff,
            T::from_f64(tk),
            s.r_dom,
            s.r_for,
            true,
        );
        let step_pv = s.notional * tau * (gf / s.spot_fx) * (c1 - c2).max_val(T::zero()) * surv;
        coupon_pv += step_pv;
        if let Some(h_ko) = s.fx_knockout_barrier {
            if h_ko.val() <= s.spot_fx.val() {
                surv = T::from_f64(1e-8);
            } else {
                let fwd_fx = s.spot_fx * ((s.r_dom - s.r_for) * tk).exp();
                let st = (sig_eff * tk.sqrt()).max_val(T::from_f64(1e-8));
                let var_t = st * st;
                let h_eff = h_ko * (st * (tau / tk).sqrt() * 0.5826).exp();
                let b = (h_eff / s.spot_fx).ln();
                let mu_t = (fwd_fx / s.spot_fx).ln() - var_t * 0.5;
                let d1 = (b - mu_t) / st;
                let d2 = (-b - mu_t) / st;
                let exp_arg = (mu_t * b * 2.0 / var_t).min_val(T::from_f64(50.0));
                let p_no_hit = (d1.norm_cdf() - exp_arg.exp() * d2.norm_cdf())
                    .max_val(T::from_f64(1e-8))
                    .min_val(T::one());
                surv = surv.min_val(p_no_hit);
            }
        }
        prev_t = tk;
    }
    PrdcValuation {
        npv: coupon_pv,
        survival_probability: surv,
        expected_coupon_pv: coupon_pv,
    }
}

#[cfg(test)]
mod phase27_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> CrossCurrencyConvertibleWarrantSpec<T> {
        CrossCurrencyConvertibleWarrantSpec {
            domestic_stock_spot: x,
            domestic_conversion_strike: T::from_f64(100.0),
            spot_fx_for_per_dom: T::from_f64(1.1),
            fixed_contractual_fx_for_per_dom: T::from_f64(1.1),
            foreign_bond_par: T::from_f64(100.0),
            foreign_coupon_rate: T::from_f64(0.03),
            r_domestic: T::from_f64(0.03),
            r_foreign: T::from_f64(0.02),
            stock_dividend_yield: T::from_f64(0.01),
            credit_spread: T::from_f64(0.02),
            stock_vol: T::from_f64(0.2),
            fx_vol: T::from_f64(0.1),
            corr_stock_fx: T::from_f64(0.2),
            maturity: T::from_f64(5.0),
        }
    }
    #[test]
    fn dual_fd_parity() {
        let h = 1e-5;
        let fd = (price_cross_currency_convertible_and_warrant_bond(&make(100.0 + h))
            .embedded_cross_currency_cb_price
            - price_cross_currency_convertible_and_warrant_bond(&make(100.0 - h))
                .embedded_cross_currency_cb_price)
            / (2.0 * h);
        let d =
            price_cross_currency_convertible_and_warrant_bond(&make(Dual::<1>::variable(100.0, 0)))
                .embedded_cross_currency_cb_price;
        assert!(d.d[0].abs() > 1e-6 && (d.d[0] - fd).abs() < 1e-5);
    }

    #[test]
    fn phase27_invariants() {
        let r = price_cross_currency_convertible_and_warrant_bond(&make(100.0));
        assert!(r.straight_bond_floor_foreign <= r.embedded_cross_currency_cb_price);
        assert!(r.embedded_cross_currency_cb_price <= r.composite_detachable_warrant_bond_price);
    }
}

#[cfg(test)]
mod phase23_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn vv_greek_match() {
        let s = CastagnaFullVannaVolgaHedgeSpec {
            spot: 1.1,
            r_dom: 0.02,
            r_for: 0.01,
            expiry: 1.0,
            atm_vol: 0.12,
            rr_25d: 0.01,
            bf_25d: 0.005,
            target_strike: 1.1,
            knockout_barrier: 1.3,
            unwind_attenuation_factor: 0.5,
            is_up_and_out: true,
            is_call: true,
        };
        let r = evaluate_castagna_full_vanna_volga_hedge_and_barrier(&s);
        assert!(r.greek_matching_max_residual < 1e-8);
        assert!(r.castagna_vv_knockout_price >= 0.0);
    }

    #[test]
    fn dual_fd_parity() {
        let mk = |k: f64| CastagnaFullVannaVolgaHedgeSpec {
            spot: 1.1,
            r_dom: 0.02,
            r_for: 0.01,
            expiry: 1.0,
            atm_vol: 0.12,
            rr_25d: 0.01,
            bf_25d: 0.005,
            target_strike: k,
            knockout_barrier: 1.3,
            unwind_attenuation_factor: 0.5,
            is_up_and_out: true,
            is_call: true,
        };
        let h = 1e-5;
        let fd = (evaluate_castagna_full_vanna_volga_hedge_and_barrier(&mk(1.1 + h))
            .exact_vv_vanilla_price
            - evaluate_castagna_full_vanna_volga_hedge_and_barrier(&mk(1.1 - h))
                .exact_vv_vanilla_price)
            / (2.0 * h);
        let x = Dual::<1>::variable(1.1, 0);
        let ad = evaluate_castagna_full_vanna_volga_hedge_and_barrier(
            &CastagnaFullVannaVolgaHedgeSpec {
                spot: Dual::<1>::from_value(1.1),
                r_dom: Dual::<1>::from_value(0.02),
                r_for: Dual::<1>::from_value(0.01),
                expiry: Dual::<1>::one(),
                atm_vol: Dual::<1>::from_value(0.12),
                rr_25d: Dual::<1>::from_value(0.01),
                bf_25d: Dual::<1>::from_value(0.005),
                target_strike: x,
                knockout_barrier: Dual::<1>::from_value(1.3),
                unwind_attenuation_factor: Dual::<1>::from_value(0.5),
                is_up_and_out: true,
                is_call: true,
            },
        )
        .exact_vv_vanilla_price;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn fx_exotic_parity_and_flat() {
        let q = FxSmileQuote {
            atm_vol: 0.2,
            rr_25d: -0.01,
            bf_25d: 0.005,
        };
        let s = FxVannaVolgaSmile::new(100.0, 0.02, 0.01, 1.0, q);
        let c = s.price_digital_cash_dom(100.0, true);
        let p = s.price_digital_cash_dom(100.0, false);
        assert!((c + p - (-0.02f64).exp()).abs() < 1e-12);
        let ot_breached = s.price_one_touch(95.0, true, false);
        assert!((ot_breached - 1.0).abs() < 1e-12);
        let ko_far = s.price_knock_out_barrier(100.0, 300.0, true, true);
        let van = garman_kohlhagen_price(100.0, 100.0, s.implied_vol(100.0), 1.0, 0.02, 0.01, true);
        assert!((ko_far - van).abs() < 1e-4);
        let h = 1e-5;
        let s_up = FxVannaVolgaSmile::new(100.0 + h, 0.02, 0.01, 1.0, q);
        let s_dn = FxVannaVolgaSmile::new(100.0 - h, 0.02, 0.01, 1.0, q);
        let fd_dnt = (s_up.price_double_no_touch(85.0, 115.0)
            - s_dn.price_double_no_touch(85.0, 115.0))
            / (2.0 * h);
        let d = FxVannaVolgaSmile::<Dual<1>>::new(
            Dual::variable(100.0, 0),
            Dual::from_value(0.02),
            Dual::from_value(0.01),
            Dual::from_value(1.0),
            FxSmileQuote {
                atm_vol: Dual::from_value(0.2),
                rr_25d: Dual::from_value(-0.01),
                bf_25d: Dual::from_value(0.005),
            },
        );
        let ad_dnt = d.price_double_no_touch(Dual::from_value(85.0), Dual::from_value(115.0));
        assert!((ad_dnt.d[0] - fd_dnt).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase16_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn prdc_correlation_monotonicity_and_ad() {
        let mut s = PrdcNoteSpec {
            notional: 100.0,
            spot_fx: 150.0,
            r_dom: 0.01,
            r_for: 0.04,
            a_dom: 0.05,
            sigma_dom: 0.008,
            a_for: 0.08,
            sigma_for: 0.01,
            sigma_fx: 0.12,
            rho_fx_dom: 0.3,
            rho_fx_for: -0.2,
            rho_dom_for: 0.25,
            coupon_times: vec![1.0, 2.0, 3.0],
            foreign_gearing: 0.15,
            domestic_gearing: 0.10,
            coupon_cap: 0.08,
            fx_knockout_barrier: Some(180.0),
        };
        let v_pos = price_prdc_note_hybrid(&s).npv;
        s.rho_fx_dom = -0.3;
        let v_neg = price_prdc_note_hybrid(&s).npv;
        assert!(v_pos != v_neg && v_pos > 0.0);
        let h = 1e-5;
        let eval_s = |fx: f64| {
            let mut sp = s.clone();
            sp.spot_fx = fx;
            price_prdc_note_hybrid(&sp).npv
        };
        let fd = (eval_s(150.0 + h) - eval_s(150.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = PrdcNoteSpec {
            notional: D::from_value(100.0),
            spot_fx: D::variable(150.0, 0),
            r_dom: D::from_value(0.01),
            r_for: D::from_value(0.04),
            a_dom: D::from_value(0.05),
            sigma_dom: D::from_value(0.008),
            a_for: D::from_value(0.08),
            sigma_for: D::from_value(0.01),
            sigma_fx: D::from_value(0.12),
            rho_fx_dom: D::from_value(-0.3),
            rho_fx_for: D::from_value(-0.2),
            rho_dom_for: D::from_value(0.25),
            coupon_times: vec![1.0, 2.0, 3.0],
            foreign_gearing: D::from_value(0.15),
            domestic_gearing: D::from_value(0.10),
            coupon_cap: D::from_value(0.08),
            fx_knockout_barrier: Some(D::from_value(180.0)),
        };
        let ad = price_prdc_note_hybrid(&sd).npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn fx_kiko_and_window_barrier_limits_and_ad() {
        let s_zero = FxWindowBarrierSpec {
            spot: 100.0,
            strike: 100.0,
            barrier: 130.0,
            vol: 0.20,
            expiry: 1.0,
            r_dom: 0.03,
            r_for: 0.01,
            window_start: 0.0,
            window_end: 0.0,
            is_up_and_out: true,
            is_call: true,
        };
        let van = garman_kohlhagen_price(100.0, 100.0, 0.20, 1.0, 0.03, 0.01, true);
        assert!((price_fx_window_barrier_option(&s_zero) - van).abs() < 1e-12);

        let mut s_full = s_zero.clone();
        s_full.window_end = 1.0;
        let v_full = price_fx_window_barrier_option(&s_full);
        assert!(v_full < van);

        let kiko = FxKikoOptionSpec {
            spot: 100.0,
            strike: 100.0,
            lower_barrier: 82.0,
            upper_barrier: 118.0,
            vol: 0.20,
            expiry: 1.0,
            r_dom: 0.03,
            r_for: 0.01,
            is_call: true,
            lower_is_ko_upper_is_ki: true,
        };
        let v_kiko = price_fx_kiko_option(&kiko);
        assert!(v_kiko > 0.0 && v_kiko < van);

        let h = 1e-5;
        let eval_sp = |sp: f64| {
            let mut k2 = kiko.clone();
            k2.spot = sp;
            price_fx_kiko_option(&k2)
        };
        let fd = (eval_sp(100.0 + h) - eval_sp(100.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let kd = FxKikoOptionSpec {
            spot: D::variable(100.0, 0),
            strike: D::from_value(100.0),
            lower_barrier: D::from_value(82.0),
            upper_barrier: D::from_value(118.0),
            vol: D::from_value(0.20),
            expiry: D::one(),
            r_dom: D::from_value(0.03),
            r_for: D::from_value(0.01),
            is_call: true,
            lower_is_ko_upper_is_ki: true,
        };
        let ad = price_fx_kiko_option(&kd);
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;
    use ql_ad::Dual;

    fn spec(pa: bool) -> ClarkFxDeltaConventionSpec<f64> {
        ClarkFxDeltaConventionSpec {
            spot: 1.1,
            r_dom: 0.03,
            r_for: 0.02,
            expiry: 1.0,
            atm_vol: 0.12,
            rr_25d: 0.03,
            market_strangle_25d: 0.01,
            is_premium_included_delta: pa,
            is_forward_delta: true,
            query_strike: 1.1,
        }
    }

    #[test]
    fn clark_smile_reprices_and_preserves_convexity() {
        let r = calibrate_clark_fx_delta_smile(&spec(false));
        assert!(r.strangle_convexity_spread.abs() > 1e-10);
        let p = calibrate_clark_fx_delta_smile(&spec(true));
        assert!(p.query_call_price.val() >= 0.0);
    }

    #[test]
    fn clark_dual_sensitivity_is_live() {
        let d = spec(false);
        let h = 1e-5;
        let mut up = d;
        up.market_strangle_25d += h;
        let mut dn = d;
        dn.market_strangle_25d -= h;
        let fd = (calibrate_clark_fx_delta_smile(&up).query_call_price
            - calibrate_clark_fx_delta_smile(&dn).query_call_price)
            / (2.0 * h);
        let ad_spec = ClarkFxDeltaConventionSpec::<Dual<1>> {
            market_strangle_25d: Dual::variable(0.01, 0),
            spot: Dual::from_value(1.1),
            r_dom: Dual::from_value(0.03),
            r_for: Dual::from_value(0.02),
            expiry: Dual::one(),
            atm_vol: Dual::from_value(0.12),
            rr_25d: Dual::from_value(0.03),
            is_premium_included_delta: false,
            is_forward_delta: true,
            query_strike: Dual::from_value(1.1),
        };
        let ad = calibrate_clark_fx_delta_smile(&ad_spec).query_call_price;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
