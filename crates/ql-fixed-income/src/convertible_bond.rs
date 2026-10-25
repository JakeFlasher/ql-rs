#![allow(clippy::needless_range_loop)]
use ql_ad::{Dual, Real};
#[derive(Clone, Debug)]
pub struct ConvertibleBondSpec<T: Real> {
    pub spot: T,
    pub conversion_ratio: T,
    pub face_value: T,
    pub coupon_rate: T,
    pub coupon_times: Vec<T>,
    pub maturity: T,
    pub rate: T,
    pub div_yield: T,
    pub credit_spread: T,
    pub vol: T,
    pub call_price: Option<T>,
    pub put_price: Option<T>,
}
pub fn price_convertible_bond_tf<T: Real>(s: &ConvertibleBondSpec<T>, steps: usize) -> (T, T) {
    let n = steps.max(1);
    let dt = s.maturity / n as f64;
    // In the zero-volatility limit the CRR up/down factors coincide and the
    // usual risk-neutral probability is 0/0.  Collapse to a deterministic
    // forward growth factor instead; using `p=0.5` keeps both branches equal
    // while preserving the correct drift.
    let deterministic = s.vol.val().abs() < 1e-14;
    let u = if deterministic {
        ((s.rate - s.div_yield) * dt).exp()
    } else {
        (s.vol * dt.sqrt()).exp()
    };
    let d = if deterministic { u } else { T::one() / u };
    let p = if deterministic {
        T::from_f64(0.5)
    } else {
        (((s.rate - s.div_yield) * dt).exp() - d) / (u - d)
    };
    let mut v = vec![T::zero(); n + 1];
    let mut b = v.clone();
    let mut step_coupons = vec![T::zero(); n + 1];
    let mut previous = T::zero();
    for &ct in &s.coupon_times {
        if ct.val() > previous.val() {
            let step_idx = ((ct / dt).val().round() as usize).clamp(1, n);
            step_coupons[step_idx] += s.face_value * s.coupon_rate * (ct - previous);
            previous = ct;
        }
    }
    let maturity_coupon = step_coupons[n];
    for i in 0..=n {
        let spot = s.spot * u.powf((n - i) as f64) * d.powf(i as f64);
        let conv = s.conversion_ratio * spot;
        let cash = s.face_value + maturity_coupon;
        v[i] = conv.max_val(cash);
        b[i] = if conv.val() >= cash.val() {
            T::zero()
        } else {
            cash
        };
    }
    for k in (0..n).rev() {
        for i in 0..=k {
            let ev = p * v[i] + (T::one() - p) * v[i + 1];
            let eb = p * b[i] + (T::one() - p) * b[i + 1];
            let mut val =
                (-s.rate * dt).exp() * (ev - eb) + (-(s.rate + s.credit_spread) * dt).exp() * eb;
            let mut bb = (-(s.rate + s.credit_spread) * dt).exp() * eb;
            if let Some(c) = s.call_price
                && val.val() > c.val()
            {
                val = c;
                bb = c;
            }
            if let Some(q) = s.put_price
                && val.val() < q.val()
            {
                val = q;
                bb = q;
            }
            if k > 0 {
                let coupon = step_coupons[k];
                val += coupon;
                bb += coupon;
            }
            let spot = s.spot * u.powf((k - i) as f64) * d.powf(i as f64);
            let conv = s.conversion_ratio * spot;
            if conv.val() > val.val() {
                val = conv;
                bb = T::zero();
            }
            v[i] = val;
            b[i] = bb;
        }
    }
    (v[0], b[0])
}
pub fn convertible_bond_greeks_dual4(
    s: &ConvertibleBondSpec<f64>,
    steps: usize,
) -> (f64, [f64; 4]) {
    let d: ConvertibleBondSpec<Dual<4>> = ConvertibleBondSpec {
        spot: Dual::<4>::variable(s.spot, 0),
        conversion_ratio: Dual::from_value(s.conversion_ratio),
        face_value: Dual::from_value(s.face_value),
        coupon_rate: Dual::from_value(s.coupon_rate),
        coupon_times: s
            .coupon_times
            .iter()
            .map(|&x| Dual::from_value(x))
            .collect(),
        maturity: Dual::from_value(s.maturity),
        rate: Dual::variable(s.rate, 2),
        div_yield: Dual::from_value(s.div_yield),
        credit_spread: Dual::variable(s.credit_spread, 3),
        vol: Dual::variable(s.vol, 1),
        call_price: s.call_price.map(Dual::from_value),
        put_price: s.put_price.map(Dual::from_value),
    };
    let (v, _) = price_convertible_bond_tf(&d, steps);
    (v.val, [v.d[0], v.d[1], v.d[2], v.d[3]])
}

#[derive(Clone, Debug)]
pub struct AfvConvertibleBondSpec<T: Real = f64> {
    pub spot: T,
    pub conversion_ratio: T,
    pub face_value: T,
    pub coupon_rate: T,
    pub maturity: T,
    pub rate: T,
    pub div_yield: T,
    pub vol: T,
    pub base_hazard_rate: T,
    pub hazard_equity_elasticity: T,
    pub equity_loss_on_default_eta: T,
    pub bond_recovery_rate: T,
    pub call_price: Option<T>,
    pub soft_call_trigger_spot: Option<T>,
    pub put_price: Option<T>,
}

#[derive(Copy, Clone, Debug)]
pub struct AfvConvertibleValuation<T: Real> {
    pub total_price: T,
    pub bond_floor: T,
    pub equity_option_value: T,
    pub local_hazard_rate_spot: T,
}

fn solve_tridiag_generic<T: Real>(l: &[T], d: &[T], u: &[T], rhs: &[T]) -> Vec<T> {
    let n = d.len();
    let mut c_p = vec![T::zero(); n];
    let mut d_p = vec![T::zero(); n];
    c_p[0] = u[0] / d[0];
    d_p[0] = rhs[0] / d[0];
    for i in 1..n {
        let denom = d[i] - l[i] * c_p[i - 1];
        if i + 1 < n {
            c_p[i] = u[i] / denom;
        }
        d_p[i] = (rhs[i] - l[i] * d_p[i - 1]) / denom;
    }
    let mut x = vec![T::zero(); n];
    x[n - 1] = d_p[n - 1];
    for i in (0..(n - 1)).rev() {
        x[i] = d_p[i] - c_p[i] * x[i + 1];
    }
    x
}

pub fn price_convertible_bond_afv_pde<T: Real>(
    spec: &AfvConvertibleBondSpec<T>,
    n_spot_nodes: usize,
    n_time_steps: usize,
) -> AfvConvertibleValuation<T> {
    let m = (n_spot_nodes.max(15)) | 1;
    let mid = m / 2;
    let nt = n_time_steps.max(8);
    let dt = spec.maturity.max_val(T::from_f64(1e-6)) / nt as f64;
    let x_max = 1.6_f64;
    let dx = 2.0 * x_max / (m - 1) as f64;
    let eps = T::from_f64(1e-4);
    let smooth_max = |a: T, b: T| -> T {
        let d = a - b;
        (a + b + (d * d + eps * eps).sqrt()) * 0.5
    };
    let smooth_min = |a: T, b: T| -> T {
        let d = a - b;
        (a + b - (d * d + eps * eps).sqrt()) * 0.5
    };

    let mut spots = Vec::with_capacity(m);
    let mut lambdas = Vec::with_capacity(m);
    for i in 0..m {
        let x_i = -x_max + i as f64 * dx;
        let s_i = spec.spot * x_i.exp();
        let lam_i = spec.base_hazard_rate * (-spec.hazard_equity_elasticity * x_i).exp();
        spots.push(s_i);
        lambdas.push(lam_i);
    }

    let mut u_cb = vec![T::zero(); m];
    let mut u_bf = vec![T::zero(); m];
    let coupon_dt = spec.face_value * spec.coupon_rate * dt;
    for i in 0..m {
        let conv = spec.conversion_ratio * spots[i];
        let redemption = spec.face_value;
        u_cb[i] = smooth_max(redemption, conv);
        u_bf[i] = redemption;
    }

    let sig2 = spec.vol * spec.vol;
    for _step in 0..nt {
        let mut l = vec![T::zero(); m];
        let mut d = vec![T::zero(); m];
        let mut u = vec![T::zero(); m];
        let mut rhs_cb = vec![T::zero(); m];
        let mut rhs_bf = vec![T::zero(); m];

        for i in 0..m {
            let lam = lambdas[i];
            let s_post_def = spots[i] * (T::one() - spec.equity_loss_on_default_eta);
            let rec_debt = spec.bond_recovery_rate * spec.face_value;
            let def_payoff_cb = smooth_max(rec_debt, spec.conversion_ratio * s_post_def);

            if i == 0 || i + 1 == m {
                d[i] = T::one() + (spec.rate + lam) * dt;
                rhs_cb[i] = u_cb[i] + coupon_dt + lam * def_payoff_cb * dt;
                rhs_bf[i] = u_bf[i] + coupon_dt + lam * rec_debt * dt;
            } else {
                // AFV (2003) Eq. 12: risk-neutral drift includes +\lambda(S) \eta for default jump martingale compensation
                let mu_x =
                    spec.rate - spec.div_yield + lam * spec.equity_loss_on_default_eta - sig2 * 0.5;
                let a_diff = sig2 * 0.5 / (dx * dx);
                let b_adv = mu_x / (2.0 * dx);
                l[i] = -(a_diff - b_adv) * dt;
                u[i] = -(a_diff + b_adv) * dt;
                d[i] = T::one() + (spec.rate + lam + a_diff * 2.0) * dt;
                rhs_cb[i] = u_cb[i] + coupon_dt + lam * def_payoff_cb * dt;
                rhs_bf[i] = u_bf[i] + coupon_dt + lam * rec_debt * dt;
            }
        }

        u_cb = solve_tridiag_generic(&l, &d, &u, &rhs_cb);
        u_bf = solve_tridiag_generic(&l, &d, &u, &rhs_bf);

        for i in 0..m {
            let conv = spec.conversion_ratio * spots[i];
            let mut val = u_cb[i];
            if let Some(call_p) = spec.call_price {
                let trig_ok = spec
                    .soft_call_trigger_spot
                    .map(|tr| spots[i].val() >= tr.val())
                    .unwrap_or(true);
                if trig_ok {
                    let call_cap = smooth_max(call_p, conv);
                    val = smooth_min(val, call_cap);
                }
            }
            if let Some(put_p) = spec.put_price {
                val = smooth_max(val, put_p);
            }
            u_cb[i] = smooth_max(val, conv);
        }
    }

    let total_price = u_cb[mid];
    let bond_floor = u_bf[mid];
    AfvConvertibleValuation {
        total_price,
        bond_floor,
        equity_option_value: (total_price - bond_floor).max_val(T::zero()),
        local_hazard_rate_spot: lambdas[mid],
    }
}

#[cfg(test)]
mod phase19_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn afv_convertible_bond_pde_bounds_hazard_monotonicity_and_ad() {
        let spec = AfvConvertibleBondSpec {
            spot: 100.0,
            conversion_ratio: 1.0,
            face_value: 100.0,
            coupon_rate: 0.03,
            maturity: 2.0,
            rate: 0.04,
            div_yield: 0.01,
            vol: 0.25,
            base_hazard_rate: 0.02,
            hazard_equity_elasticity: 1.2,
            equity_loss_on_default_eta: 0.5,
            bond_recovery_rate: 0.40,
            call_price: Some(115.0),
            soft_call_trigger_spot: Some(130.0),
            put_price: Some(95.0),
        };
        let v = price_convertible_bond_afv_pde(&spec, 21, 16);
        assert!(v.total_price >= v.bond_floor);
        assert!(v.total_price >= spec.conversion_ratio * spec.spot - 1e-2);
        assert!(v.equity_option_value > 0.0);

        let mut high_haz = spec.clone();
        high_haz.base_hazard_rate = 0.06;
        let v_high = price_convertible_bond_afv_pde(&high_haz, 21, 16);
        assert!(v_high.bond_floor < v.bond_floor);

        let h = 1e-5;
        let eval_s = |s0: f64| {
            let mut sp = spec.clone();
            sp.spot = s0;
            price_convertible_bond_afv_pde(&sp, 21, 16).total_price
        };
        let fd = (eval_s(100.0 + h) - eval_s(100.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = AfvConvertibleBondSpec {
            spot: D::variable(100.0, 0),
            conversion_ratio: D::one(),
            face_value: D::from_value(100.0),
            coupon_rate: D::from_value(0.03),
            maturity: D::from_value(2.0),
            rate: D::from_value(0.04),
            div_yield: D::from_value(0.01),
            vol: D::from_value(0.25),
            base_hazard_rate: D::from_value(0.02),
            hazard_equity_elasticity: D::from_value(1.2),
            equity_loss_on_default_eta: D::from_value(0.5),
            bond_recovery_rate: D::from_value(0.40),
            call_price: Some(D::from_value(115.0)),
            soft_call_trigger_spot: Some(D::from_value(130.0)),
            put_price: Some(D::from_value(95.0)),
        };
        let ad = price_convertible_bond_afv_pde(&spec_d, 21, 16).total_price;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-4);
    }
}

#[derive(Clone, Debug)]
pub struct SirbuShreveDilutionMakeWholeConvertibleSpec<T: Real = f64> {
    pub spot: T,
    pub volatility: T,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub credit_spread: T,
    pub maturity_years: T,
    pub par_value: T,
    pub coupon_rate_annual: T,
    pub conversion_ratio: T,
    pub existing_shares_millions: T,
    pub convertible_issue_count_millions: T,
    pub debt_relief_retention_fraction: T,
    pub call_price: T,
    pub soft_call_trigger_spot: T,
    pub make_whole_max_bonus_shares: T,
    pub make_whole_time_decay: T,
    pub make_whole_moneyness_decay: T,
    pub tree_steps: usize,
}
#[derive(Clone, Debug)]
pub struct SirbuShreveDilutionMakeWholeConvertibleReport<T: Real = f64> {
    pub dilution_factor_q: T,
    pub diluted_convertible_price: T,
    pub undiluted_convertible_price: T,
    pub dilution_price_discount: T,
    pub make_whole_enhanced_price: T,
    pub make_whole_option_value: T,
    pub effective_conversion_ratio_t0: T,
    pub diluted_parity_t0: T,
    pub equity_delta: T,
}
fn smooth_pos<T: Real>(x: T) -> T {
    (x + (x * x + T::from_f64(1e-8)).sqrt()) / T::from_f64(2.0)
}
fn smooth_max<T: Real>(a: T, b: T) -> T {
    a + smooth_pos(b - a)
}

fn mw_bonus_shares<T: Real>(
    s: &SirbuShreveDilutionMakeWholeConvertibleSpec<T>,
    sp: T,
    t_frac: T,
) -> T {
    let s_soft = s.soft_call_trigger_spot.max_val(T::from_f64(1e-12));
    let excess = smooth_pos(sp - s.soft_call_trigger_spot) / s_soft;
    s.make_whole_max_bonus_shares
        * (-s.make_whole_time_decay * t_frac).exp()
        * (-s.make_whole_moneyness_decay * excess).exp()
}

fn cb_tree<T: Real>(s: &SirbuShreveDilutionMakeWholeConvertibleSpec<T>, q: T, mw: bool) -> (T, T) {
    let n = s.tree_steps.max(1);
    let dt = s.maturity_years / T::from_f64(n as f64);
    let u = (s.volatility * dt.sqrt()).exp();
    let d = T::one() / u;
    let den = (u - d).max_val(T::from_f64(1e-12));
    let p = (((s.risk_free_rate - s.dividend_yield) * dt).exp() - d) / den;
    let disc = (-(s.risk_free_rate + s.credit_spread * T::from_f64(0.5)) * dt).exp();
    let coup = s.par_value * s.coupon_rate_annual * dt;
    let mut v = vec![T::zero(); n + 1];
    for i in 0..=n {
        let sp = s.spot * u.powf((n - i) as f64) * d.powf(i as f64);
        let g = s.conversion_ratio
            + if mw {
                mw_bonus_shares(s, sp, T::one())
            } else {
                T::zero()
            };
        let parity = (T::one() - q) * g * sp + q * s.par_value * s.debt_relief_retention_fraction;
        v[i] = parity.max_val(s.par_value);
    }
    let mut delta_t0 = T::zero();
    for k in (0..n).rev() {
        let mut next = vec![T::zero(); k + 1];
        let t_frac = T::from_f64(k as f64 / n as f64);
        for i in 0..=k {
            let sp = s.spot * u.powf((k - i) as f64) * d.powf(i as f64);
            let ev = disc * (p * v[i] + (T::one() - p) * v[i + 1] + coup);
            let g = s.conversion_ratio
                + if mw {
                    mw_bonus_shares(s, sp, t_frac)
                } else {
                    T::zero()
                };
            let parity =
                (T::one() - q) * g * sp + q * s.par_value * s.debt_relief_retention_fraction;
            let w = T::one()
                / (T::one()
                    + (-(sp - s.soft_call_trigger_spot)
                        / (T::from_f64(0.05)
                            * s.soft_call_trigger_spot.max_val(T::from_f64(1e-12))))
                    .exp());
            let call = smooth_max(s.call_price, parity);
            next[i] = (ev * (T::one() - w) + call * w).max_val(parity);
        }
        if k == 1 {
            let sp_up = s.spot * u;
            let sp_dn = s.spot * d;
            delta_t0 = (next[0] - next[1]) / (s.spot * (u - d)).max_val(T::from_f64(1e-12));
            let _ = (sp_up, sp_dn);
        }
        v = next;
    }
    (v[0], delta_t0)
}

pub fn price_sirbu_shreve_dilution_and_makewhole_convertible<T: Real>(
    s: &SirbuShreveDilutionMakeWholeConvertibleSpec<T>,
) -> SirbuShreveDilutionMakeWholeConvertibleReport<T> {
    let q = (s.convertible_issue_count_millions * s.conversion_ratio)
        / (s.existing_shares_millions + s.convertible_issue_count_millions * s.conversion_ratio)
            .max_val(T::from_f64(1e-12));
    let (und, _) = cb_tree(s, T::zero(), false);
    let (dil, _) = cb_tree(s, q, false);
    let (enh_raw, delta) = cb_tree(s, q, true);
    let enh = enh_raw.max_val(dil);
    let er = s.conversion_ratio + mw_bonus_shares(s, s.spot, T::zero());
    SirbuShreveDilutionMakeWholeConvertibleReport {
        dilution_factor_q: q,
        diluted_convertible_price: dil,
        undiluted_convertible_price: und,
        dilution_price_discount: und - dil,
        make_whole_enhanced_price: enh,
        make_whole_option_value: (enh - dil).max_val(T::zero()),
        effective_conversion_ratio_t0: er,
        diluted_parity_t0: (T::one() - q) * er * s.spot
            + q * s.par_value * s.debt_relief_retention_fraction,
        equity_delta: delta,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ThorpKassoufCalamosConvertibleArbSpec<T: Real = f64> {
    pub stock_price: T,
    pub conversion_ratio: T,
    pub convertible_market_price: T,
    pub face_value: T,
    pub annual_coupon_rate: T,
    pub stock_dividend_yield: T,
    pub risk_free_rate: T,
    pub credit_spread: T,
    pub short_stock_rebate_rate: T,
    pub implied_volatility: T,
    pub realized_volatility: T,
    pub maturity: T,
}
#[derive(Clone, Copy, Debug)]
pub struct ThorpKassoufCalamosConvertibleArbReport<T: Real = f64> {
    pub conversion_parity: T,
    pub conversion_premium_ratio: T,
    pub straight_bond_investment_floor: T,
    pub theoretical_convertible_value: T,
    pub delta_neutral_short_shares: T,
    pub convertible_gamma: T,
    pub annual_cash_flow_advantage: T,
    pub breakeven_payback_years: T,
    pub expected_annual_gamma_carry_pnl: T,
}
pub fn analyze_thorp_kassouf_calamos_convertible_arbitrage<T: Real>(
    s: &ThorpKassoufCalamosConvertibleArbSpec<T>,
) -> ThorpKassoufCalamosConvertibleArbReport<T> {
    let t = s.maturity.max_val(T::from_f64(1e-12));
    let parity = s.stock_price * s.conversion_ratio;
    let y = (s.risk_free_rate + s.credit_spread).max_val(T::from_f64(1e-8));
    let coupon_pv = s.annual_coupon_rate * s.face_value * (T::one() - (-y * t).exp()) / y;
    let floor = coupon_pv + s.face_value * (-y * t).exp();
    let vol = s.implied_volatility.max_val(T::from_f64(1e-8));
    let d1 = ((s.stock_price / (s.face_value / s.conversion_ratio + T::from_f64(1e-14))).ln()
        + (s.risk_free_rate - s.stock_dividend_yield + T::from_f64(0.5) * vol * vol) * t)
        / (vol * t.sqrt());
    let d2 = d1 - vol * t.sqrt();
    let call = s.stock_price * (-s.stock_dividend_yield * t).exp() * d1.norm_cdf()
        - (s.face_value / s.conversion_ratio) * (-(s.risk_free_rate) * t).exp() * d2.norm_cdf();
    let theo = floor + call * s.conversion_ratio;
    let delta = s.conversion_ratio * (-s.stock_dividend_yield * t).exp() * d1.norm_cdf();
    let gamma = s.conversion_ratio * (-s.stock_dividend_yield * t).exp() * d1.norm_pdf()
        / (s.stock_price * vol * t.sqrt());
    let adv = s.annual_coupon_rate * s.face_value
        + s.short_stock_rebate_rate * s.stock_price * delta
        - s.stock_dividend_yield * s.stock_price * delta
        - s.credit_spread * floor;
    let pay = (s.convertible_market_price - parity).abs() / (adv.abs() + T::from_f64(1e-14));
    let carry = gamma
        * s.stock_price
        * s.stock_price
        * (s.realized_volatility * s.realized_volatility
            - s.implied_volatility * s.implied_volatility)
        / T::from_f64(2.0)
        + adv;
    ThorpKassoufCalamosConvertibleArbReport {
        conversion_parity: parity,
        conversion_premium_ratio: s.convertible_market_price / parity.max_val(T::from_f64(1e-14))
            - T::one(),
        straight_bond_investment_floor: floor,
        theoretical_convertible_value: theo,
        delta_neutral_short_shares: delta,
        convertible_gamma: gamma,
        annual_cash_flow_advantage: adv,
        breakeven_payback_years: pay,
        expected_annual_gamma_carry_pnl: carry,
    }
}

#[cfg(test)]
mod phase24_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: Real>(x: T) -> SirbuShreveDilutionMakeWholeConvertibleSpec<T> {
        SirbuShreveDilutionMakeWholeConvertibleSpec {
            spot: x,
            volatility: T::from_f64(0.2),
            risk_free_rate: T::from_f64(0.03),
            dividend_yield: T::zero(),
            credit_spread: T::from_f64(0.02),
            maturity_years: T::from_f64(3.),
            par_value: T::from_f64(100.),
            coupon_rate_annual: T::from_f64(0.03),
            conversion_ratio: T::from_f64(0.8),
            existing_shares_millions: T::from_f64(100.),
            convertible_issue_count_millions: T::from_f64(10.),
            debt_relief_retention_fraction: T::zero(),
            call_price: T::from_f64(105.),
            soft_call_trigger_spot: T::from_f64(110.),
            make_whole_max_bonus_shares: T::from_f64(0.2),
            make_whole_time_decay: T::from_f64(0.2),
            make_whole_moneyness_decay: T::from_f64(0.3),
            tree_steps: 16,
        }
    }
    #[test]
    fn phase24() {
        let r = price_sirbu_shreve_dilution_and_makewhole_convertible(&mk(100.));
        assert!(r.undiluted_convertible_price > r.diluted_convertible_price);
        assert!(r.make_whole_enhanced_price >= r.diluted_convertible_price);
        let h = 1e-5;
        let fd = (price_sirbu_shreve_dilution_and_makewhole_convertible(&mk(100. + h))
            .make_whole_enhanced_price
            - price_sirbu_shreve_dilution_and_makewhole_convertible(&mk(100. - h))
                .make_whole_enhanced_price)
            / (2. * h);
        let ad = price_sirbu_shreve_dilution_and_makewhole_convertible(&mk(Dual::<1>::variable(
            100., 0,
        )))
        .make_whole_enhanced_price;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase25_tests {
    use super::*;
    use ql_ad::Dual;
    fn cb<T: Real>(x: T) -> ThorpKassoufCalamosConvertibleArbSpec<T> {
        ThorpKassoufCalamosConvertibleArbSpec {
            stock_price: x,
            conversion_ratio: T::from_f64(0.8),
            convertible_market_price: T::from_f64(95.0),
            face_value: T::from_f64(100.0),
            annual_coupon_rate: T::from_f64(0.03),
            stock_dividend_yield: T::from_f64(0.01),
            risk_free_rate: T::from_f64(0.02),
            credit_spread: T::from_f64(0.03),
            short_stock_rebate_rate: T::from_f64(0.01),
            implied_volatility: T::from_f64(0.2),
            realized_volatility: T::from_f64(0.22),
            maturity: T::from_f64(3.0),
        }
    }
    #[test]
    fn convertible_arb_dual_fd() {
        let r = analyze_thorp_kassouf_calamos_convertible_arbitrage(&cb(100.0));
        assert!(
            (r.conversion_parity * (1.0 + r.conversion_premium_ratio)
                - cb(100.0).convertible_market_price)
                .abs()
                < 1e-12
        );
        let h = 1e-5;
        let fd = (analyze_thorp_kassouf_calamos_convertible_arbitrage(&cb(100.0 + h))
            .theoretical_convertible_value
            .val()
            - analyze_thorp_kassouf_calamos_convertible_arbitrage(&cb(100.0 - h))
                .theoretical_convertible_value
                .val())
            / (2.0 * h);
        let ad =
            analyze_thorp_kassouf_calamos_convertible_arbitrage(&cb(Dual::<1>::variable(100.0, 0)))
                .theoretical_convertible_value;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
