#![allow(clippy::needless_range_loop)]
use ql_ad::{Dual, Real};
use ql_vol::black_scholes_price;
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

#[derive(Clone, Copy, Debug)]
pub struct KoziolAsquithFasbStrategicConvertibleSpec<T: Real = f64> {
    pub firm_value: T,
    pub debt_face_value: T,
    pub existing_shares: T,
    pub new_shares_blockholder1: T,
    pub new_shares_blockholder2: T,
    pub coupon_rate: T,
    pub corporate_tax_rate: T,
    pub dividend_yield: T,
    pub stock_price: T,
    pub base_conversion_ratio: T,
    pub induced_sweetened_conversion_ratio: T,
    pub cash_sweetener_amount: T,
    pub call_price: T,
    pub firm_volatility: T,
    pub risk_free_rate: T,
    pub maturity: T,
}
#[derive(Clone, Copy, Debug)]
pub struct KoziolAsquithFasbStrategicConvertibleReport<T: Real = f64> {
    pub blockholder1_dilution_factor: T,
    pub joint_conversion_dilution_factor: T,
    pub koziol_preemptive_threshold_firm_value: T,
    pub koziol_cooperative_threshold_firm_value: T,
    pub prisoners_dilemma_early_exercise_wedge: T,
    pub asquith_after_tax_call_delay_advantage: T,
    pub asquith_breakeven_dividend_call_trigger: T,
    pub fasb_asu2024_04_induced_conversion_expense: T,
    pub strategic_convertible_fair_value: T,
}
pub fn evaluate_koziol_asquith_fasb_strategic_convertible<T: Real>(
    s: &KoziolAsquithFasbStrategicConvertibleSpec<T>,
) -> KoziolAsquithFasbStrategicConvertibleReport<T> {
    let total = (s.existing_shares + s.new_shares_blockholder1 + s.new_shares_blockholder2)
        .max_val(T::from_f64(1e-10));
    let alpha1_solo = s.new_shares_blockholder1
        / (s.existing_shares + s.new_shares_blockholder1).max_val(T::from_f64(1e-10));
    let alpha1_joint = s.new_shares_blockholder1 / total;
    let alpha2_joint = s.new_shares_blockholder2 / total;
    let d1 = T::one() - alpha1_solo;
    let dj = T::one() - alpha1_joint - alpha2_joint;
    let sig2 = s.firm_volatility * s.firm_volatility;
    let beta_base =
        T::from_f64(0.5) - (s.risk_free_rate - s.dividend_yield) / sig2.max_val(T::from_f64(1e-10));
    let beta = beta_base
        + (beta_base * beta_base
            + T::from_f64(2.0) * s.risk_free_rate / sig2.max_val(T::from_f64(1e-10)))
        .sqrt();
    let coop =
        beta / (beta - T::one()).max_val(T::from_f64(1e-10)) * T::from_f64(2.0) * s.debt_face_value
            / (alpha1_joint + alpha2_joint).max_val(T::from_f64(1e-10));
    let pre = coop
        * (alpha1_joint / alpha1_solo.max_val(T::from_f64(1e-10)))
            .max_val(T::from_f64(1e-10))
            .pow(T::one() / beta.max_val(T::from_f64(1e-10)));
    let wedge = (coop - pre).abs();
    let annuity = (T::one() - (-s.risk_free_rate * s.maturity).exp())
        / s.risk_free_rate.max_val(T::from_f64(1e-10));
    let tax_adv = (s.dividend_yield * s.base_conversion_ratio * s.stock_price
        - (T::one() - s.corporate_tax_rate) * s.coupon_rate * s.debt_face_value)
        * annuity;
    let trigger = (T::one() - s.corporate_tax_rate) * s.coupon_rate * s.debt_face_value
        / (s.base_conversion_ratio * s.stock_price).max_val(T::from_f64(1e-10));
    let expense = (s.induced_sweetened_conversion_ratio - s.base_conversion_ratio)
        .max_val(T::zero())
        * s.stock_price
        + s.cash_sweetener_amount;
    let option = black_scholes_price(
        s.stock_price,
        s.call_price,
        s.firm_volatility.max_val(T::from_f64(1e-8)),
        s.maturity.max_val(T::from_f64(1e-8)),
        s.risk_free_rate,
        s.dividend_yield,
        true,
    );
    let leverage_adj =
        (s.debt_face_value / s.firm_value.max_val(T::from_f64(1e-10))).min_val(T::one());
    let floor = (s.debt_face_value * (-s.risk_free_rate * s.maturity).exp()
        + s.coupon_rate * s.debt_face_value * annuity)
        * (T::one() - T::from_f64(0.01) * leverage_adj);
    let fair = floor + option + expense.max_val(T::zero());
    KoziolAsquithFasbStrategicConvertibleReport {
        blockholder1_dilution_factor: d1,
        joint_conversion_dilution_factor: dj,
        koziol_preemptive_threshold_firm_value: pre,
        koziol_cooperative_threshold_firm_value: coop,
        prisoners_dilemma_early_exercise_wedge: wedge,
        asquith_after_tax_call_delay_advantage: tax_adv,
        asquith_breakeven_dividend_call_trigger: trigger,
        fasb_asu2024_04_induced_conversion_expense: expense.max_val(T::zero()),
        strategic_convertible_fair_value: fair,
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
mod domain_suite_36_tests {
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
    fn domain_suite_36() {
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
mod domain_suite_37_tests {
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

#[cfg(test)]
mod domain_suite_41_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn koziol_invariants_and_dual() {
        fn mk<T: Real>(x: T) -> KoziolAsquithFasbStrategicConvertibleSpec<T> {
            KoziolAsquithFasbStrategicConvertibleSpec {
                firm_value: T::from_f64(10000.0),
                debt_face_value: T::from_f64(8000.0),
                existing_shares: T::from_f64(100.0),
                new_shares_blockholder1: T::from_f64(10.0),
                new_shares_blockholder2: T::from_f64(8.0),
                coupon_rate: T::from_f64(0.04),
                corporate_tax_rate: T::from_f64(0.25),
                dividend_yield: T::from_f64(0.01),
                stock_price: x,
                base_conversion_ratio: T::one(),
                induced_sweetened_conversion_ratio: T::from_f64(1.1),
                cash_sweetener_amount: T::one(),
                call_price: T::from_f64(120.0),
                firm_volatility: T::from_f64(0.2),
                risk_free_rate: T::from_f64(0.03),
                maturity: T::from_f64(5.0),
            }
        }
        let r = evaluate_koziol_asquith_fasb_strategic_convertible(&mk(100.0));
        assert!(
            r.blockholder1_dilution_factor > r.joint_conversion_dilution_factor
                && r.koziol_preemptive_threshold_firm_value
                    < r.koziol_cooperative_threshold_firm_value
                && r.fasb_asu2024_04_induced_conversion_expense >= 0.0
        );
        let h = 1e-5;
        let fd = (evaluate_koziol_asquith_fasb_strategic_convertible(&mk(100.0 + h))
            .strategic_convertible_fair_value
            - evaluate_koziol_asquith_fasb_strategic_convertible(&mk(100.0 - h))
                .strategic_convertible_fair_value)
            / (2.0 * h);
        let ad =
            evaluate_koziol_asquith_fasb_strategic_convertible(&mk(Dual::<1>::variable(100.0, 0)))
                .strategic_convertible_fair_value;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_31_tests {
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

mod evaluate_an_daoquan_clifford_chance_resettable_convertible_module {
    use ql_ad::Real;
    #[derive(Clone, Copy, Debug)]
    pub struct AnDaoquanCliffordChanceResettableCbSpec<T: Real = f64> {
        pub spot_price: T,
        pub initial_conversion_price: T,
        pub book_value_per_share_bps: T,
        pub avg_20day_price: T,
        pub face_value: T,
        pub coupon_rate: T,
        pub redemption_maturity_price: T,
        pub put_back_trigger_ratio: T,
        pub put_back_price: T,
        pub soft_call_trigger_ratio: T,
        pub soft_call_price: T,
        pub stock_volatility: T,
        pub risk_free_rate: T,
        pub credit_spread: T,
        pub maturity: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct AnDaoquanCliffordChanceResettableCbReport<T: Real = f64> {
        pub straight_bond_floor: T,
        pub initial_conversion_parity: T,
        pub legally_floored_reset_conversion_price: T,
        pub post_reset_conversion_parity: T,
        pub downward_reset_option_uplift: T,
        pub conditional_put_back_protection_value: T,
        pub issuer_soft_call_cap_giveup: T,
        pub three_line_duplex_fair_value: T,
    }
    pub fn evaluate_an_daoquan_clifford_chance_resettable_convertible<T: Real>(
        s: &AnDaoquanCliffordChanceResettableCbSpec<T>,
    ) -> AnDaoquanCliffordChanceResettableCbReport<T> {
        let floor = (s.redemption_maturity_price
            + s.face_value * s.coupon_rate * s.maturity * T::from_f64(0.5))
            * (-(s.risk_free_rate + s.credit_spread) * s.maturity).exp();
        let init_price = s.initial_conversion_price.max_val(T::from_f64(1e-6));
        let initial = s.face_value / init_price * s.spot_price;
        let candidate = s.avg_20day_price.max_val(s.book_value_per_share_bps);
        let post_price = init_price.min_val(candidate).max_val(T::from_f64(1e-6));
        let post = s.face_value / post_price * s.spot_price;
        let trigger_gap = s.initial_conversion_price * s.put_back_trigger_ratio - s.avg_20day_price;
        let reset_weight = T::one() / (T::one() + (-trigger_gap * T::from_f64(0.2)).exp());
        let uplift = (post - initial).max_val(T::zero())
            * (T::from_f64(0.5) + T::from_f64(0.5) * reset_weight)
            + s.stock_volatility * s.maturity.sqrt() * T::from_f64(4.0);
        let put = (s.put_back_price - initial.min_val(floor)).max_val(T::zero())
            * (-s.risk_free_rate * s.maturity).exp()
            * (T::from_f64(0.25) + T::from_f64(0.5) * reset_weight);
        let soft_call_threshold = s
            .soft_call_price
            .max_val(s.initial_conversion_price * s.soft_call_trigger_ratio);
        let cap = (post - soft_call_threshold).max_val(T::zero()) * T::from_f64(0.2);
        let fair = floor.max_val(initial) + uplift + put - cap;
        AnDaoquanCliffordChanceResettableCbReport {
            straight_bond_floor: floor,
            initial_conversion_parity: initial,
            legally_floored_reset_conversion_price: post_price,
            post_reset_conversion_parity: post,
            downward_reset_option_uplift: uplift,
            conditional_put_back_protection_value: put,
            issuer_soft_call_cap_giveup: cap,
            three_line_duplex_fair_value: fair,
        }
    }

    #[cfg(test)]
    mod domain_suite_42_tests {
        use super::*;
        use ql_ad::Dual;
        fn make<T: Real>(x: T) -> AnDaoquanCliffordChanceResettableCbSpec<T> {
            AnDaoquanCliffordChanceResettableCbSpec {
                spot_price: T::from_f64(90.0) + x,
                initial_conversion_price: T::from_f64(100.0),
                book_value_per_share_bps: T::from_f64(80.0),
                avg_20day_price: T::from_f64(85.0),
                face_value: T::from_f64(100.0),
                coupon_rate: T::from_f64(0.02),
                redemption_maturity_price: T::from_f64(105.0),
                put_back_trigger_ratio: T::from_f64(0.7),
                put_back_price: T::from_f64(103.0),
                soft_call_trigger_ratio: T::from_f64(1.3),
                soft_call_price: T::from_f64(105.0),
                stock_volatility: T::from_f64(0.2),
                risk_free_rate: T::from_f64(0.03),
                credit_spread: T::from_f64(0.01),
                maturity: T::from_f64(3.0),
            }
        }
        #[test]
        fn suite42_reset_cb_dual_fd() {
            let r = evaluate_an_daoquan_clifford_chance_resettable_convertible(&make(0.0));
            assert!(r.legally_floored_reset_conversion_price >= 80.0);
            assert!(r.post_reset_conversion_parity > r.initial_conversion_parity);
            let h = 1e-6;
            let fd = (evaluate_an_daoquan_clifford_chance_resettable_convertible(&make(h))
                .post_reset_conversion_parity
                - evaluate_an_daoquan_clifford_chance_resettable_convertible(&make(-h))
                    .post_reset_conversion_parity)
                / (2.0 * h);
            let ad = evaluate_an_daoquan_clifford_chance_resettable_convertible(&make(
                Dual::<1>::variable(0.0, 0),
            ))
            .post_reset_conversion_parity;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_an_daoquan_clifford_chance_resettable_convertible_module::*;

mod evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole_module {
    use ql_ad::Real;
    #[derive(Clone, Copy, Debug)]
    pub struct MayerBrownSsgaCappedCallMakeWholeSpec<T: Real = f64> {
        pub bond_par: T,
        pub spot_price: T,
        pub base_conversion_price: T,
        pub cap_strike_price: T,
        pub volatility: T,
        pub risk_free_rate: T,
        pub maturity_years: T,
        pub acquisition_offer_price: T,
        pub remaining_years_at_acquisition: T,
        pub tax_rate_section_163j: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct MayerBrownSsgaCappedCallMakeWholeReport<T: Real = f64> {
        pub base_conversion_ratio: T,
        pub lower_call_option_value: T,
        pub upper_cap_call_option_value: T,
        pub capped_call_overlay_gross_cost: T,
        pub after_tax_integrated_capped_call_cost: T,
        pub unhedged_net_dilution_shares_at_offer: T,
        pub capped_call_hedged_dilution_shares_at_offer: T,
        pub make_whole_additional_shares: T,
        pub total_fundamental_change_conversion_value: T,
    }
    fn call<T: Real>(s: T, k: T, v: T, r: T, t: T) -> T {
        let vol = (v * t.sqrt()).max_val(T::from_f64(1e-8));
        let d1 = ((s / k).ln() + (r + T::from_f64(0.5) * v * v) * t) / vol;
        let d2 = d1 - vol;
        s * d1.norm_cdf() - k * (-r * t).exp() * d2.norm_cdf()
    }
    pub fn evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole<T: Real>(
        s: &MayerBrownSsgaCappedCallMakeWholeSpec<T>,
    ) -> MayerBrownSsgaCappedCallMakeWholeReport<T> {
        let ratio = s.bond_par / s.base_conversion_price.max_val(T::from_f64(1e-8));
        let lower = ratio
            * call(
                s.spot_price,
                s.base_conversion_price,
                s.volatility,
                s.risk_free_rate,
                s.maturity_years,
            );
        let upper = ratio
            * call(
                s.spot_price,
                s.cap_strike_price,
                s.volatility,
                s.risk_free_rate,
                s.maturity_years,
            );
        let gross = (lower - upper).max_val(T::zero());
        let after = gross * (T::one() - s.tax_rate_section_163j);
        let s_offer_safe = s.acquisition_offer_price.max_val(T::from_f64(1e-8));
        let unhedged = ratio
            * ((s.acquisition_offer_price - s.base_conversion_price).max_val(T::zero())
                / s_offer_safe);
        let hedged = ratio
            * ((s.acquisition_offer_price - s.cap_strike_price).max_val(T::zero()) / s_offer_safe);
        let intrinsic =
            (s.acquisition_offer_price - s.base_conversion_price).max_val(T::zero()) * ratio;
        let call_at_offer = ratio
            * call(
                s.acquisition_offer_price,
                s.base_conversion_price,
                s.volatility,
                s.risk_free_rate,
                s.remaining_years_at_acquisition,
            );
        let mw = (call_at_offer - intrinsic).max_val(T::zero()) / s_offer_safe;
        MayerBrownSsgaCappedCallMakeWholeReport {
            base_conversion_ratio: ratio,
            lower_call_option_value: lower,
            upper_cap_call_option_value: upper,
            capped_call_overlay_gross_cost: gross,
            after_tax_integrated_capped_call_cost: after,
            unhedged_net_dilution_shares_at_offer: unhedged,
            capped_call_hedged_dilution_shares_at_offer: hedged,
            make_whole_additional_shares: mw,
            total_fundamental_change_conversion_value: (ratio + mw) * s.acquisition_offer_price,
        }
    }
    #[cfg(test)]
    mod domain_suite_44_tests {
        use super::*;
        use ql_ad::Dual;

        fn spec<T: Real>(x: T) -> MayerBrownSsgaCappedCallMakeWholeSpec<T> {
            MayerBrownSsgaCappedCallMakeWholeSpec {
                bond_par: T::from_f64(100.0),
                spot_price: T::from_f64(100.0) + x,
                base_conversion_price: T::from_f64(90.0),
                cap_strike_price: T::from_f64(120.0),
                volatility: T::from_f64(0.2),
                risk_free_rate: T::from_f64(0.03),
                maturity_years: T::from_f64(5.0),
                acquisition_offer_price: T::from_f64(110.0),
                remaining_years_at_acquisition: T::from_f64(2.0),
                tax_rate_section_163j: T::from_f64(0.2),
            }
        }

        #[test]
        fn overlay_positive() {
            let r = evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(&spec(0.0));
            assert!(r.capped_call_overlay_gross_cost > 0.0);
            assert!(
                r.capped_call_hedged_dilution_shares_at_offer
                    <= r.unhedged_net_dilution_shares_at_offer
            );
            assert!(r.make_whole_additional_shares > 0.0);
            let h = 1e-6;
            let fd = (evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(&spec(h))
                .capped_call_overlay_gross_cost
                - evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(&spec(-h))
                    .capped_call_overlay_gross_cost)
                / (2.0 * h);
            let ad = evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole(&spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .capped_call_overlay_gross_cost;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_mayer_brown_ssga_capped_call_overlay_and_make_whole_module::*;

mod evaluate_despiegeleer_mandatory_convertible_decs_and_percs_module {
    use ql_ad::Real;
    #[derive(Clone, Copy, Debug)]
    pub struct DeSpiegeleerMandatoryConvertibleDecsSpec<T: Real = f64> {
        pub spot_price: T,
        pub par_amount: T,
        pub lower_strike_kl: T,
        pub upper_strike_ku: T,
        pub decs_coupon_rate: T,
        pub risk_free_rate: T,
        pub dividend_yield: T,
        pub credit_spread: T,
        pub lower_strike_vol: T,
        pub upper_strike_vol: T,
        pub maturity_years: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct DeSpiegeleerMandatoryConvertibleDecsReport<T: Real = f64> {
        pub max_conversion_ratio_lower: T,
        pub min_conversion_ratio_upper: T,
        pub fixed_coupon_strip_pv: T,
        pub prepaid_forward_equity_pv: T,
        pub short_lower_call_leg_pv: T,
        pub long_upper_call_leg_pv: T,
        pub percs_capped_fair_value: T,
        pub decs_mandatory_fair_value: T,
        pub decs_over_percs_upside_participation: T,
        pub decomposition_conservation_error: T,
    }
    fn call<T: Real>(s: T, k: T, v: T, r: T, q: T, t: T) -> T {
        let st = t.sqrt();
        let d1 = ((s / k).ln() + (r - q + v * v / T::from_f64(2.0)) * t) / (v * st);
        let d2 = d1 - v * st;
        s * (-q * t).exp() * d1.norm_cdf() - k * (-r * t).exp() * d2.norm_cdf()
    }
    pub fn evaluate_despiegeleer_mandatory_convertible_decs_and_percs<T: Real>(
        s: &DeSpiegeleerMandatoryConvertibleDecsSpec<T>,
    ) -> DeSpiegeleerMandatoryConvertibleDecsReport<T> {
        let nmax = s.par_amount / s.lower_strike_kl;
        let nmin = s.par_amount / s.upper_strike_ku;
        let z = s.risk_free_rate + s.credit_spread;
        let coup = s.decs_coupon_rate * s.par_amount * (T::one() - (-z * s.maturity_years).exp())
            / z.max_val(T::from_f64(1e-12));
        let fwd = nmax * s.spot_price * (-s.dividend_yield * s.maturity_years).exp();
        let lo = nmax
            * call(
                s.spot_price,
                s.lower_strike_kl,
                s.lower_strike_vol,
                s.risk_free_rate,
                s.dividend_yield,
                s.maturity_years,
            );
        let hi = nmin
            * call(
                s.spot_price,
                s.upper_strike_ku,
                s.upper_strike_vol,
                s.risk_free_rate,
                s.dividend_yield,
                s.maturity_years,
            );
        let p = coup + fwd - lo;
        let d = p + hi;
        DeSpiegeleerMandatoryConvertibleDecsReport {
            max_conversion_ratio_lower: nmax,
            min_conversion_ratio_upper: nmin,
            fixed_coupon_strip_pv: coup,
            prepaid_forward_equity_pv: fwd,
            short_lower_call_leg_pv: -lo,
            long_upper_call_leg_pv: hi,
            percs_capped_fair_value: p,
            decs_mandatory_fair_value: d,
            decs_over_percs_upside_participation: hi,
            decomposition_conservation_error: (d - (coup + fwd - lo + hi)).abs(),
        }
    }
    #[derive(Clone, Copy, Debug)]
    pub struct DavidsonLevinFourFactorIoPoSpec<T: Real = f64> {
        pub initial_balance: T,
        pub gross_wac_rate: T,
        pub net_passthrough_coupon: T,
        pub current_mortgage_rate: T,
        pub base_turnover_cpr: T,
        pub max_refi_cpr: T,
        pub burnout_sensitivity_beta: T,
        pub seasonality_amplitude: T,
        pub oas_discount_rate: T,
        pub num_periods: usize,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct DavidsonLevinFourFactorIoPoReport<T: Real = f64> {
        pub terminal_pool_factor: T,
        pub average_realized_cpr: T,
        pub passthrough_present_value: T,
        pub io_strip_present_value: T,
        pub po_strip_present_value: T,
        pub io_po_strip_conservation_error: T,
        pub io_effective_rate_sensitivity: T,
        pub po_effective_rate_sensitivity: T,
    }
    pub fn evaluate_davidson_levin_four_factor_prepayment_io_po<T: Real>(
        s: &DavidsonLevinFourFactorIoPoSpec<T>,
    ) -> DavidsonLevinFourFactorIoPoReport<T> {
        let mut bal = T::one();
        let mut dbal_dr = T::zero();
        let mut io = T::zero();
        let mut po = T::zero();
        let mut dio_dr = T::zero();
        let mut dpo_dr = T::zero();
        let n_per = s.num_periods.max(1);
        for i in 1..=n_per {
            let t_yr = T::from_f64(i as f64 / 12.0);
            let age = T::from_f64((i as f64 / 30.0).min(1.0));
            let bal_safe = bal.max_val(T::from_f64(1e-12));
            let burnout = (bal_safe.ln() * s.burnout_sensitivity_beta).exp();
            let dburn_dbal = s.burnout_sensitivity_beta * burnout / bal_safe;
            let exp_arg = (T::from_f64(-40.0)
                * (s.gross_wac_rate - s.current_mortgage_rate - T::from_f64(0.005)))
            .exp();
            let sig = T::one() / (T::one() + exp_arg);
            let dsig_drmtg = -T::from_f64(40.0) * exp_arg * sig * sig;
            let refi = s.max_refi_cpr * burnout * sig;
            let drefi_dr = s.max_refi_cpr * (dburn_dbal * dbal_dr * sig + burnout * dsig_drmtg);
            let seas = T::one()
                + s.seasonality_amplitude
                    * T::from_f64((i as f64 * std::f64::consts::PI / 6.0).sin());
            let cpr_annual = ((s.base_turnover_cpr * age + refi) * seas)
                .min_val(T::from_f64(0.95))
                .max_val(T::from_f64(1e-6));
            let dcpr_dr = drefi_dr * seas;
            let one_minus_cpr = (T::one() - cpr_annual).max_val(T::from_f64(1e-6));
            let smm = T::one() - one_minus_cpr.powf(1.0 / 12.0);
            let dsmm_dr = (T::one() - smm) / (T::from_f64(12.0) * one_minus_cpr) * dcpr_dr;
            let disc = (-s.oas_discount_rate * t_yr).exp();
            let c_m = s.net_passthrough_coupon / T::from_f64(12.0);
            io += bal * c_m * disc;
            dio_dr += dbal_dr * c_m * disc;
            let pre = if i == n_per { bal } else { bal * smm };
            let dpre_dr = if i == n_per {
                dbal_dr
            } else {
                dbal_dr * smm + bal * dsmm_dr
            };
            po += pre * disc;
            dpo_dr += dpre_dr * disc;
            let next_bal = bal * (T::one() - smm);
            dbal_dr = dbal_dr * (T::one() - smm) - bal * dsmm_dr;
            bal = next_bal;
        }
        let pv = io + po;
        let avg = T::one()
            - (bal.max_val(T::from_f64(1e-12)).ln() * T::from_f64(12.0 / (n_per as f64))).exp();
        let err = (pv - io - po) * s.initial_balance;
        DavidsonLevinFourFactorIoPoReport {
            terminal_pool_factor: bal,
            average_realized_cpr: avg,
            passthrough_present_value: pv * s.initial_balance,
            io_strip_present_value: io * s.initial_balance,
            po_strip_present_value: po * s.initial_balance,
            io_po_strip_conservation_error: (err * err + T::from_f64(1e-30)).sqrt(),
            io_effective_rate_sensitivity: dio_dr * s.initial_balance,
            po_effective_rate_sensitivity: dpo_dr * s.initial_balance,
        }
    }

    #[cfg(test)]
    mod domain_suite_46_tests {
        use super::*;
        use ql_ad::Dual;
        fn decs_spec<T: Real>(x: T) -> DeSpiegeleerMandatoryConvertibleDecsSpec<T> {
            DeSpiegeleerMandatoryConvertibleDecsSpec {
                spot_price: T::from_f64(100.0) + x,
                par_amount: T::from_f64(100.0),
                lower_strike_kl: T::from_f64(90.0),
                upper_strike_ku: T::from_f64(120.0),
                decs_coupon_rate: T::from_f64(0.06),
                risk_free_rate: T::from_f64(0.03),
                dividend_yield: T::from_f64(0.02),
                credit_spread: T::from_f64(0.015),
                lower_strike_vol: T::from_f64(0.25),
                upper_strike_vol: T::from_f64(0.22),
                maturity_years: T::from_f64(3.0),
            }
        }
        fn mbs_spec<T: Real>(x: T) -> DavidsonLevinFourFactorIoPoSpec<T> {
            DavidsonLevinFourFactorIoPoSpec {
                initial_balance: T::from_f64(100.0),
                gross_wac_rate: T::from_f64(0.06),
                net_passthrough_coupon: T::from_f64(0.055),
                current_mortgage_rate: T::from_f64(0.05) + x,
                base_turnover_cpr: T::from_f64(0.06),
                max_refi_cpr: T::from_f64(0.35),
                burnout_sensitivity_beta: T::from_f64(0.5),
                seasonality_amplitude: T::from_f64(0.1),
                oas_discount_rate: T::from_f64(0.045),
                num_periods: 12,
            }
        }
        #[test]
        fn suite46_fixed_income_invariants_and_dual_fd() {
            let h = 1e-6;
            let d0 = evaluate_despiegeleer_mandatory_convertible_decs_and_percs(&decs_spec(0.0));
            assert!(
                d0.decomposition_conservation_error < 1e-12
                    && d0.decs_mandatory_fair_value > d0.percs_capped_fair_value
            );
            let fd_d = (evaluate_despiegeleer_mandatory_convertible_decs_and_percs(&decs_spec(h))
                .decs_mandatory_fair_value
                - evaluate_despiegeleer_mandatory_convertible_decs_and_percs(&decs_spec(-h))
                    .decs_mandatory_fair_value)
                / (2.0 * h);
            let ad_d = evaluate_despiegeleer_mandatory_convertible_decs_and_percs(&decs_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .decs_mandatory_fair_value;
            assert!(ad_d.d[0].abs() > 1e-6 && (ad_d.d[0] - fd_d).abs() < 1e-5);

            let m0 = evaluate_davidson_levin_four_factor_prepayment_io_po(&mbs_spec(0.0));
            assert!(
                m0.io_po_strip_conservation_error < 1e-12
                    && m0.io_effective_rate_sensitivity > 0.0
                    && m0.po_effective_rate_sensitivity < 0.0
            );
            let fd_m = (evaluate_davidson_levin_four_factor_prepayment_io_po(&mbs_spec(h))
                .passthrough_present_value
                - evaluate_davidson_levin_four_factor_prepayment_io_po(&mbs_spec(-h))
                    .passthrough_present_value)
                / (2.0 * h);
            let ad_m = evaluate_davidson_levin_four_factor_prepayment_io_po(&mbs_spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .passthrough_present_value;
            assert!(ad_m.d[0].abs() > 1e-6 && (ad_m.d[0] - fd_m).abs() < 1e-4);
        }
    }
}
pub use evaluate_despiegeleer_mandatory_convertible_decs_and_percs_module::*;

mod evaluate_despiegeleer_cb_soft_call_make_whole_module {
    // De Spiegeleer, Schoutens & Van Hulle (2014) Convertible Bond Soft-Call Parity Trigger,
    // Make-Whole Additional Shares Table & Stock-Borrow Fee Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct DeSpiegeleerCbSoftCallMakeWholeSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub par_value: T,
        pub conversion_ratio_cr: T,
        pub soft_call_trigger_multiple: T,
        pub make_whole_max_bonus_shares: T,
        pub coupon_rate_c: T,
        pub risky_bond_yield_y: T,
        pub risk_free_rate_r: T,
        pub dividend_and_borrow_cost_q: T,
        pub equity_volatility_sigma: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct DeSpiegeleerCbSoftCallMakeWholeReport<T: Real = f64> {
        pub conversion_price_k_conv: T,
        pub soft_call_barrier_price_b: T,
        pub conversion_parity_value: T,
        pub straight_bond_floor: T,
        pub embedded_conversion_call_value: T,
        pub soft_call_issuer_cap_reduction: T,
        pub make_whole_table_protection_pv: T,
        pub convertible_bond_fair_value: T,
    }

    fn bs_call<T: Real>(s: T, k: T, r: T, q: T, sig: T, t: T) -> T {
        let std = sig * t.sqrt();
        let d1 = ((s / k).ln() + (r - q + T::from_f64(0.5) * sig * sig) * t) / std;
        let d2 = d1 - std;
        s * (-q * t).exp() * norm_cdf(d1) - k * (-r * t).exp() * norm_cdf(d2)
    }

    pub fn evaluate_despiegeleer_cb_soft_call_make_whole<T: Real>(
        spec: &DeSpiegeleerCbSoftCallMakeWholeSpec<T>,
    ) -> DeSpiegeleerCbSoftCallMakeWholeReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));
        let par = spec.par_value.max_val(T::one());
        let cr = spec.conversion_ratio_cr.max_val(T::from_f64(1e-4));
        let trig_mult = spec.soft_call_trigger_multiple.max_val(T::from_f64(1.05));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));
        let sig = spec.equity_volatility_sigma.max_val(T::from_f64(1e-4));
        let y = spec.risky_bond_yield_y.max_val(T::from_f64(1e-4));

        let k_conv = par / cr;
        let b_soft = k_conv * trig_mult;
        let parity = cr * s0;

        // Straight risky bond floor:
        let annuity = (T::one() - (-y * t).exp()) / y;
        let bond_floor = par * spec.coupon_rate_c * annuity + par * (-y * t).exp();

        // Embedded conversion call struck at effective floor strike K_eff = BondFloor / CR:
        let k_eff = (bond_floor / cr).max_val(T::from_f64(1e-4));
        let conv_call = cr
            * bs_call(
                s0,
                k_eff,
                spec.risk_free_rate_r,
                spec.dividend_and_borrow_cost_q,
                sig,
                t,
            );

        // Provisional soft-call cap above B_soft (weighted by M-of-N excursion damping):
        let soft_cap_raw = cr
            * bs_call(
                s0,
                b_soft,
                spec.risk_free_rate_r,
                spec.dividend_and_borrow_cost_q,
                sig,
                t,
            );
        let excursion_prob = (-T::from_f64(0.35) * sig * t.sqrt()).exp();
        let soft_call_reduction = T::from_f64(0.5) * excursion_prob * soft_cap_raw;

        // Make-whole bonus shares protection PV when called or acquired above K_conv:
        let mw_pv = spec.make_whole_max_bonus_shares
            * (bs_call(
                s0,
                k_conv,
                spec.risk_free_rate_r,
                spec.dividend_and_borrow_cost_q,
                sig,
                t,
            ) - bs_call(
                s0,
                b_soft,
                spec.risk_free_rate_r,
                spec.dividend_and_borrow_cost_q,
                sig,
                t,
            ))
            .max_val(T::zero());

        let cb_val = bond_floor + conv_call - soft_call_reduction + mw_pv;

        DeSpiegeleerCbSoftCallMakeWholeReport {
            conversion_price_k_conv: k_conv,
            soft_call_barrier_price_b: b_soft,
            conversion_parity_value: parity,
            straight_bond_floor: bond_floor,
            embedded_conversion_call_value: conv_call,
            soft_call_issuer_cap_reduction: soft_call_reduction,
            make_whole_table_protection_pv: mw_pv,
            convertible_bond_fair_value: cb_val,
        }
    }

    #[cfg(test)]
    mod domain_suite_54_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn despiegeleer_cb_soft_call_make_whole_invariants_and_ad_parity() {
            let spec = DeSpiegeleerCbSoftCallMakeWholeSpec {
                spot_price_s0: 105.0,
                par_value: 100.0,
                conversion_ratio_cr: 1.0,
                soft_call_trigger_multiple: 1.30,
                make_whole_max_bonus_shares: 0.18,
                coupon_rate_c: 0.025,
                risky_bond_yield_y: 0.055,
                risk_free_rate_r: 0.040,
                dividend_and_borrow_cost_q: 0.018,
                equity_volatility_sigma: 0.28,
                maturity_years_t: 4.0,
            };
            let rep = evaluate_despiegeleer_cb_soft_call_make_whole(&spec);
            assert!(rep.convertible_bond_fair_value > rep.straight_bond_floor);
            assert!(rep.soft_call_barrier_price_b > rep.conversion_price_k_conv);
            assert!(rep.make_whole_table_protection_pv > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_despiegeleer_cb_soft_call_make_whole(
                &DeSpiegeleerCbSoftCallMakeWholeSpec {
                    spot_price_s0: Dual::<1>::variable(spec.spot_price_s0, 0),
                    par_value: Dual::<1>::from_value(spec.par_value),
                    conversion_ratio_cr: Dual::<1>::from_value(spec.conversion_ratio_cr),
                    soft_call_trigger_multiple: Dual::<1>::from_value(
                        spec.soft_call_trigger_multiple,
                    ),
                    make_whole_max_bonus_shares: Dual::<1>::from_value(
                        spec.make_whole_max_bonus_shares,
                    ),
                    coupon_rate_c: Dual::<1>::from_value(spec.coupon_rate_c),
                    risky_bond_yield_y: Dual::<1>::from_value(spec.risky_bond_yield_y),
                    risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
                    dividend_and_borrow_cost_q: Dual::<1>::from_value(
                        spec.dividend_and_borrow_cost_q,
                    ),
                    equity_volatility_sigma: Dual::<1>::from_value(spec.equity_volatility_sigma),
                    maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                },
            );
            let mut sp = spec;
            sp.spot_price_s0 += h;
            let mut sm = spec;
            sm.spot_price_s0 -= h;
            let fd = (evaluate_despiegeleer_cb_soft_call_make_whole(&sp)
                .convertible_bond_fair_value
                - evaluate_despiegeleer_cb_soft_call_make_whole(&sm).convertible_bond_fair_value)
                / (2.0 * h);
            assert!((r_ad.convertible_bond_fair_value.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_despiegeleer_cb_soft_call_make_whole_module::*;

mod evaluate_tsiveriotis_fernandes_split_discount_convertible_module {
    // Tsiveriotis & Fernandes (1998, J. Fixed Income) & De Spiegeleer-Schoutens (2011, Ch. 7)
    // Credit-Blended Convertible Bond (Split Cash-Only Debt vs. Equity-Conversion Discounting) Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct TsiveriotisFernandesSplitConvertibleSpec<T: Real = f64> {
        pub stock_price_s0: T,
        pub conversion_ratio_cr: T,
        pub par_redemption_value: T,
        pub annual_coupon_rate_c: T,
        pub risk_free_rate_r: T,
        pub credit_spread_sc: T,
        pub dividend_yield_q: T,
        pub equity_volatility_sigma: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct TsiveriotisFernandesSplitConvertibleReport<T: Real = f64> {
        pub parity_conversion_value: T,
        pub risky_straight_bond_floor: T,
        pub risk_neutral_conversion_probability: T,
        pub equity_settled_component_pv: T,
        pub cash_only_debt_component_pv: T,
        pub tsiveriotis_fernandes_convertible_value: T,
        pub blended_effective_discount_rate: T,
        pub single_rate_risky_discount_wedge: T,
    }

    pub fn evaluate_tsiveriotis_fernandes_split_discount_convertible<T: Real>(
        spec: &TsiveriotisFernandesSplitConvertibleSpec<T>,
    ) -> TsiveriotisFernandesSplitConvertibleReport<T> {
        let s0 = spec.stock_price_s0.max_val(T::from_f64(1.0));
        let cr = spec.conversion_ratio_cr.max_val(T::from_f64(0.01));
        let par = spec.par_redemption_value.max_val(T::from_f64(10.0));
        let c = spec.annual_coupon_rate_c.max_val(T::zero());
        let r = spec.risk_free_rate_r;
        let sc = spec.credit_spread_sc.max_val(T::from_f64(1e-4));
        let q = spec.dividend_yield_q;
        let sig = spec.equity_volatility_sigma.max_val(T::from_f64(0.05));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        let parity = cr * s0;
        let r_risky = r + sc;
        let coupon_pv_risky =
            par * c * (T::one() - (-r_risky * t).exp()) / r_risky.max_val(T::from_f64(1e-4));
        let bond_floor = coupon_pv_risky + par * (-r_risky * t).exp();

        let k_conv = par / cr;
        let std = sig * t.sqrt();
        let d1 = ((s0 / k_conv).ln() + (r - q + T::from_f64(0.5) * sig * sig) * t) / std;
        let d2 = d1 - std;
        let p_conv = norm_cdf(d2);

        // Tsiveriotis-Fernandes (1998) split:
        // Equity-settled component E_0 is discounted at risk-free rate r (no credit spread on own stock delivery);
        // Cash-only debt component U_0 (non-conversion principal + coupons) is discounted at r + s_c:
        let eq_comp = cr * s0 * (-q * t).exp() * norm_cdf(d1);
        let cash_comp = coupon_pv_risky + par * (-r_risky * t).exp() * (T::one() - p_conv);
        let tf_val = eq_comp + cash_comp;

        // Naive single-rate risky discounting applies (r + s_c) to the embedded call as well:
        let call_risky = (-r_risky * t).exp()
            * (cr * s0 * ((r - q) * t).exp() * norm_cdf(d1) - par * norm_cdf(d2));
        let single_risky_val = bond_floor + call_risky;
        let wedge = tf_val - single_risky_val;

        let w_eq = eq_comp / tf_val.max_val(T::from_f64(1.0));
        let r_blended = (T::one() - w_eq) * r_risky + w_eq * r;

        TsiveriotisFernandesSplitConvertibleReport {
            parity_conversion_value: parity,
            risky_straight_bond_floor: bond_floor,
            risk_neutral_conversion_probability: p_conv,
            equity_settled_component_pv: eq_comp,
            cash_only_debt_component_pv: cash_comp,
            tsiveriotis_fernandes_convertible_value: tf_val,
            blended_effective_discount_rate: r_blended,
            single_rate_risky_discount_wedge: wedge,
        }
    }

    #[cfg(test)]
    mod evaluate_tsiveriotis_fernandes_split_discount_convertible_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn tsiveriotis_fernandes_split_convertible_invariants_and_ad_parity() {
            let spec = TsiveriotisFernandesSplitConvertibleSpec {
                stock_price_s0: 96.0,
                conversion_ratio_cr: 1.0,
                par_redemption_value: 100.0,
                annual_coupon_rate_c: 0.035,
                risk_free_rate_r: 0.040,
                credit_spread_sc: 0.025,
                dividend_yield_q: 0.012,
                equity_volatility_sigma: 0.28,
                maturity_years_t: 3.0,
            };
            let rep = evaluate_tsiveriotis_fernandes_split_discount_convertible(&spec);
            assert!(rep.tsiveriotis_fernandes_convertible_value > rep.risky_straight_bond_floor);
            assert!(rep.tsiveriotis_fernandes_convertible_value > rep.parity_conversion_value);
            assert!(rep.single_rate_risky_discount_wedge > 0.0);

            let h = 1e-5;
            let r_ad = evaluate_tsiveriotis_fernandes_split_discount_convertible(
                &TsiveriotisFernandesSplitConvertibleSpec {
                    stock_price_s0: Dual::<9>::variable(spec.stock_price_s0, 0),
                    conversion_ratio_cr: Dual::<9>::variable(spec.conversion_ratio_cr, 1),
                    par_redemption_value: Dual::<9>::variable(spec.par_redemption_value, 2),
                    annual_coupon_rate_c: Dual::<9>::variable(spec.annual_coupon_rate_c, 3),
                    risk_free_rate_r: Dual::<9>::variable(spec.risk_free_rate_r, 4),
                    credit_spread_sc: Dual::<9>::variable(spec.credit_spread_sc, 5),
                    dividend_yield_q: Dual::<9>::variable(spec.dividend_yield_q, 6),
                    equity_volatility_sigma: Dual::<9>::variable(spec.equity_volatility_sigma, 7),
                    maturity_years_t: Dual::<9>::variable(spec.maturity_years_t, 8),
                },
            );
            for (idx, &ad) in r_ad
                .tsiveriotis_fernandes_convertible_value
                .d
                .iter()
                .enumerate()
            {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.stock_price_s0 += h;
                        sm.stock_price_s0 -= h;
                    }
                    1 => {
                        sp.conversion_ratio_cr += h;
                        sm.conversion_ratio_cr -= h;
                    }
                    2 => {
                        sp.par_redemption_value += h;
                        sm.par_redemption_value -= h;
                    }
                    3 => {
                        sp.annual_coupon_rate_c += h;
                        sm.annual_coupon_rate_c -= h;
                    }
                    4 => {
                        sp.risk_free_rate_r += h;
                        sm.risk_free_rate_r -= h;
                    }
                    5 => {
                        sp.credit_spread_sc += h;
                        sm.credit_spread_sc -= h;
                    }
                    6 => {
                        sp.dividend_yield_q += h;
                        sm.dividend_yield_q -= h;
                    }
                    7 => {
                        sp.equity_volatility_sigma += h;
                        sm.equity_volatility_sigma -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_tsiveriotis_fernandes_split_discount_convertible(&sp)
                    .tsiveriotis_fernandes_convertible_value
                    - evaluate_tsiveriotis_fernandes_split_discount_convertible(&sm)
                        .tsiveriotis_fernandes_convertible_value)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_tsiveriotis_fernandes_split_discount_convertible_module::*;

mod evaluate_brennan_schwartz_callable_puttable_convertible_game_module {
    // Brennan & Schwartz (1977, Journal of Finance 32(5), "Convertible Bonds: Valuation and
    // Optimal Strategies for Call and Conversion"; 1980), Ingersoll (1977, JFE 4(3)),
    // & De Spiegeleer & Schoutens (2011, The Handbook of Convertible Bonds, Ch. 4-6)
    // Callable-Puttable Convertible Bond Game-Option Decomposition Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct BrennanSchwartzCallablePuttableCbSpec<T: Real = f64> {
        pub stock_price_s0: T,
        pub conversion_ratio_n_c: T,
        pub par_notional_n: T,
        pub coupon_rate_c: T,
        pub risk_free_rate_r: T,
        pub credit_spread_s_c: T,
        pub dividend_yield_q: T,
        pub equity_volatility_sigma_s: T,
        pub investor_put_price_k_put: T,
        pub issuer_call_trigger_price_k_call: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BrennanSchwartzCallablePuttableCbReport<T: Real = f64> {
        pub straight_risky_bond_floor: T,
        pub conversion_parity_value: T,
        pub embedded_conversion_call_value: T,
        pub investor_put_floor_option_value: T,
        pub issuer_call_forcing_giveup_value: T,
        pub game_convertible_fair_value: T,
        pub effective_equity_dollar_delta: T,
    }

    pub fn evaluate_brennan_schwartz_callable_puttable_convertible_game<T: Real>(
        spec: &BrennanSchwartzCallablePuttableCbSpec<T>,
    ) -> BrennanSchwartzCallablePuttableCbReport<T> {
        let s0 = spec.stock_price_s0.max_val(T::from_f64(1.0));
        let n_c = spec.conversion_ratio_n_c.max_val(T::from_f64(0.01));
        let par = spec.par_notional_n.max_val(T::from_f64(1.0));
        let c = spec.coupon_rate_c;
        let r = spec.risk_free_rate_r;
        let s_c = spec.credit_spread_s_c.max_val(T::from_f64(0.0005));
        let q = spec.dividend_yield_q;
        let sig_s = spec.equity_volatility_sigma_s.max_val(T::from_f64(0.02));
        let k_put = spec.investor_put_price_k_put.max_val(T::from_f64(1.0));
        let k_call = spec
            .issuer_call_trigger_price_k_call
            .max_val(par + T::from_f64(1.0));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        let y_risky = (r + s_c).max_val(T::from_f64(1e-4));
        let df_risky = (-y_risky * t).exp();
        let b_floor = par * ((c / y_risky) * (T::one() - df_risky) + df_risky);
        let parity = n_c * s0;

        let sqrt_t = t.sqrt();
        let std_s = sig_s * sqrt_t;
        let df_r = (-r * t).exp();
        let df_q = (-q * t).exp();

        // Conversion call option on n_c shares above effective per-share bond floor strike:
        let k_conv = (b_floor / n_c).max_val(T::from_f64(0.1));
        let d1_conv = ((s0 / k_conv).ln() + (r - q + T::from_f64(0.5) * sig_s * sig_s) * t) / std_s;
        let d2_conv = d1_conv - std_s;
        let c_conv = n_c * (s0 * df_q * norm_cdf(d1_conv) - k_conv * df_r * norm_cdf(d2_conv));

        // Issuer forced-conversion call option above per-share call trigger S_call* = K_call / n_c:
        let s_call_star = (k_call / n_c).max_val(k_conv + T::from_f64(0.5));
        let d1_call =
            ((s0 / s_call_star).ln() + (r - q + T::from_f64(0.5) * sig_s * sig_s) * t) / std_s;
        let d2_call = d1_call - std_s;
        // Call protection attenuates issuer call-forcing by 0.5 over the life of the bond:
        let call_giveup = T::from_f64(0.5)
            * n_c
            * (s0 * df_q * norm_cdf(d1_call) - s_call_star * df_r * norm_cdf(d2_call));

        // Investor put floor option at intermediate put date T_p = 0.5 * T on the risky bond floor:
        let t_p = T::from_f64(0.5) * t;
        let sig_b = T::from_f64(0.06) + T::from_f64(0.5) * t * s_c;
        let std_b = sig_b * t_p.sqrt();
        let d1_put = ((b_floor / k_put).ln() + T::from_f64(0.5) * sig_b * sig_b * t_p) / std_b;
        let d2_put = d1_put - std_b;
        let p_floor = (-r * t_p).exp() * (k_put * norm_cdf(-d2_put) - b_floor * norm_cdf(-d1_put));

        let cb_fair = b_floor + c_conv + p_floor - call_giveup;
        let eq_delta = n_c * df_q * (norm_cdf(d1_conv) - T::from_f64(0.5) * norm_cdf(d1_call));

        BrennanSchwartzCallablePuttableCbReport {
            straight_risky_bond_floor: b_floor,
            conversion_parity_value: parity,
            embedded_conversion_call_value: c_conv,
            investor_put_floor_option_value: p_floor,
            issuer_call_forcing_giveup_value: call_giveup,
            game_convertible_fair_value: cb_fair,
            effective_equity_dollar_delta: eq_delta,
        }
    }

    #[cfg(test)]
    mod evaluate_brennan_schwartz_callable_puttable_convertible_game_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brennan_schwartz_callable_puttable_cb_invariants_and_ad_parity() {
            let spec = BrennanSchwartzCallablePuttableCbSpec {
                stock_price_s0: 98.0,
                conversion_ratio_n_c: 1.0,
                par_notional_n: 100.0,
                coupon_rate_c: 0.025,
                risk_free_rate_r: 0.04,
                credit_spread_s_c: 0.022,
                dividend_yield_q: 0.012,
                equity_volatility_sigma_s: 0.32,
                investor_put_price_k_put: 96.0,
                issuer_call_trigger_price_k_call: 130.0,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_brennan_schwartz_callable_puttable_convertible_game(&spec);
            assert!(rep.game_convertible_fair_value > rep.straight_risky_bond_floor);
            assert!(rep.embedded_conversion_call_value > rep.issuer_call_forcing_giveup_value);
            assert!(rep.effective_equity_dollar_delta > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_brennan_schwartz_callable_puttable_convertible_game(
                &BrennanSchwartzCallablePuttableCbSpec {
                    stock_price_s0: Dual::<6>::variable(spec.stock_price_s0, 0),
                    conversion_ratio_n_c: Dual::from_f64(spec.conversion_ratio_n_c),
                    par_notional_n: Dual::from_f64(spec.par_notional_n),
                    coupon_rate_c: Dual::from_f64(spec.coupon_rate_c),
                    risk_free_rate_r: Dual::<6>::variable(spec.risk_free_rate_r, 1),
                    credit_spread_s_c: Dual::<6>::variable(spec.credit_spread_s_c, 2),
                    dividend_yield_q: Dual::from_f64(spec.dividend_yield_q),
                    equity_volatility_sigma_s: Dual::<6>::variable(
                        spec.equity_volatility_sigma_s,
                        3,
                    ),
                    investor_put_price_k_put: Dual::<6>::variable(spec.investor_put_price_k_put, 4),
                    issuer_call_trigger_price_k_call: Dual::<6>::variable(
                        spec.issuer_call_trigger_price_k_call,
                        5,
                    ),
                    maturity_years_t: Dual::from_f64(spec.maturity_years_t),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.stock_price_s0,
                    spec.risk_free_rate_r,
                    spec.credit_spread_s_c,
                    spec.equity_volatility_sigma_s,
                    spec.investor_put_price_k_put,
                    spec.issuer_call_trigger_price_k_call,
                ];
                v[i] += h;
                let up = evaluate_brennan_schwartz_callable_puttable_convertible_game(
                    &BrennanSchwartzCallablePuttableCbSpec {
                        stock_price_s0: v[0],
                        risk_free_rate_r: v[1],
                        credit_spread_s_c: v[2],
                        equity_volatility_sigma_s: v[3],
                        investor_put_price_k_put: v[4],
                        issuer_call_trigger_price_k_call: v[5],
                        ..spec
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_brennan_schwartz_callable_puttable_convertible_game(
                    &BrennanSchwartzCallablePuttableCbSpec {
                        stock_price_s0: v[0],
                        risk_free_rate_r: v[1],
                        credit_spread_s_c: v[2],
                        equity_volatility_sigma_s: v[3],
                        investor_put_price_k_put: v[4],
                        issuer_call_trigger_price_k_call: v[5],
                        ..spec
                    },
                );
                let fd =
                    (up.game_convertible_fair_value - dn.game_convertible_fair_value) / (2.0 * h);
                assert!(
                    (r_ad.game_convertible_fair_value.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_brennan_schwartz_callable_puttable_convertible_game_module::*;

/// Specification for the Tsiveriotis–Fernandes (1998) convertible bond Cash-Only vs. Equity credit split engine.
///
/// Reference: Tsiveriotis & Fernandes (1998, *Journal of Fixed Income* 8(2), "Valuing Convertible Bonds
/// with Credit Risk"), Ayache, Forsyth & Vetzal (2003), Tuckman & Serrat (2022).
#[derive(Debug, Clone, Copy)]
pub struct TsiveriotisFernandesConvertibleSpec<T: Real> {
    /// Current underlying stock price $S_0$.
    pub spot_price_s0: T,
    /// Conversion ratio $\gamma_{\text{conv}}$ (shares per bond).
    pub conversion_ratio_gamma: T,
    /// Bond redemption par value $F_{\text{par}}$.
    pub par_value: T,
    /// Annualized coupon rate $c$.
    pub coupon_rate_c: T,
    /// Risk-free interest rate $r$.
    pub risk_free_rate_r: T,
    /// Issuer credit spread $s_c$ applied to the Cash-Only Convertible Bond (COCB) component $U$.
    pub credit_spread_sc: T,
    /// Underlying stock dividend yield $q$.
    pub dividend_yield_q: T,
    /// Stock return volatility $\sigma$.
    pub equity_volatility_sigma: T,
    /// Issuer call price $K_{\text{call}}$.
    pub call_price_k: T,
    /// Maturity $T$ (in years).
    pub maturity_years_t: T,
}

/// Analytics returned by [`evaluate_tsiveriotis_fernandes_convertible_bond_credit_split`].
#[derive(Debug, Clone, Copy)]
pub struct TsiveriotisFernandesConvertibleAnalytics<T: Real> {
    /// Total convertible bond value $V_0 = U_0 + E_0$.
    pub convertible_bond_value: T,
    /// Cash-Only Convertible Bond (COCB) credit-risky component $U_0$.
    pub cash_only_cocb_component: T,
    /// Credit-risk-free equity conversion component $E_0 = V_0 - U_0$.
    pub equity_conversion_component: T,
    /// Effective blended credit discount spread $s_{\text{eff}} = s_c \cdot (U_0 / V_0)$.
    pub effective_blended_credit_spread: T,
}

/// Evaluates the Tsiveriotis–Fernandes (1998) coupled two-component $(V, U)$ backward induction on a
/// 4-step binomial lattice with smooth conversion and call exercise boundaries.
pub fn evaluate_tsiveriotis_fernandes_convertible_bond_credit_split<T: Real>(
    spec: &TsiveriotisFernandesConvertibleSpec<T>,
) -> TsiveriotisFernandesConvertibleAnalytics<T> {
    const STEPS: usize = 4;
    let t = spec.maturity_years_t.max_val(T::from_f64(1e-4));
    let dt = t / T::from_f64(STEPS as f64);
    let sig = spec.equity_volatility_sigma.max_val(T::from_f64(1e-4));
    let u_mult = (sig * dt.sqrt()).exp();
    let d_mult = T::one() / u_mult;
    let growth = ((spec.risk_free_rate_r - spec.dividend_yield_q) * dt).exp();
    let p_up = (growth - d_mult) / (u_mult - d_mult);
    let p_dn = T::one() - p_up;

    let df_rf = (-spec.risk_free_rate_r * dt).exp();
    let df_risky = (-(spec.risk_free_rate_r + spec.credit_spread_sc) * dt).exp();
    let step_coupon = spec.par_value * spec.coupon_rate_c * dt;
    let final_redemption = spec.par_value;
    let temp = spec.par_value * T::from_f64(0.06);

    let mut v_nodes = [T::zero(); STEPS + 1];
    let mut u_nodes = [T::zero(); STEPS + 1];

    // Terminal payoff at step N
    for j in 0..=STEPS {
        let up_pow = (2 * j as isize - STEPS as isize) as f64;
        let s_nj = spec.spot_price_s0 * (sig * dt.sqrt() * T::from_f64(up_pow)).exp();
        let cv_nj = spec.conversion_ratio_gamma * s_nj;
        let pi_conv = T::one() / (T::one() + (-(cv_nj - final_redemption) / temp).exp());
        v_nodes[j] = pi_conv * cv_nj + (T::one() - pi_conv) * final_redemption;
        u_nodes[j] = (T::one() - pi_conv) * final_redemption;
    }

    // Backward induction across steps n = STEPS-1 ..= 0
    for n in (0..STEPS).rev() {
        for j in 0..=n {
            let up_pow = (2 * j as isize - n as isize) as f64;
            let s_nj = spec.spot_price_s0 * (sig * dt.sqrt() * T::from_f64(up_pow)).exp();
            let cv_nj = spec.conversion_ratio_gamma * s_nj;

            // Cash-only part + period coupon discounted at risky rate r + s_c
            let u_cont = df_risky * (p_up * u_nodes[j + 1] + p_dn * u_nodes[j] + step_coupon);
            // Equity part E = V - U discounted at risk-free rate r
            let e_up = v_nodes[j + 1] - u_nodes[j + 1];
            let e_dn = v_nodes[j] - u_nodes[j];
            let e_cont = df_rf * (p_up * e_up + p_dn * e_dn);
            let v_cont = u_cont + e_cont;

            // Smooth issuer call at K_call and smooth holder conversion at cv_nj
            let pi_call = T::one() / (T::one() + (-(v_cont - spec.call_price_k) / temp).exp());
            let v_after_call = (T::one() - pi_call) * v_cont + pi_call * spec.call_price_k;
            let u_after_call = (T::one() - pi_call) * u_cont + pi_call * spec.call_price_k;

            let pi_conv = T::one() / (T::one() + (-(cv_nj - v_after_call) / temp).exp());
            v_nodes[j] = pi_conv * cv_nj + (T::one() - pi_conv) * v_after_call;
            u_nodes[j] = (T::one() - pi_conv) * u_after_call;
        }
    }

    let v0 = v_nodes[0];
    let u0 = u_nodes[0];
    let e0 = v0 - u0;
    let s_eff = spec.credit_spread_sc * (u0 / v0.max_val(T::from_f64(1e-6)));

    TsiveriotisFernandesConvertibleAnalytics {
        convertible_bond_value: v0,
        cash_only_cocb_component: u0,
        equity_conversion_component: e0,
        effective_blended_credit_spread: s_eff,
    }
}

#[cfg(test)]
mod tsiveriotis_fernandes_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn tsiveriotis_fernandes_convertible_invariants_and_ad_parity() {
        let base = TsiveriotisFernandesConvertibleSpec {
            spot_price_s0: 98.0_f64,
            conversion_ratio_gamma: 1.0,
            par_value: 100.0,
            coupon_rate_c: 0.035,
            risk_free_rate_r: 0.04,
            credit_spread_sc: 0.025,
            dividend_yield_q: 0.015,
            equity_volatility_sigma: 0.28,
            call_price_k: 112.0,
            maturity_years_t: 2.0,
        };
        let r0 = evaluate_tsiveriotis_fernandes_convertible_bond_credit_split(&base);
        assert!(r0.convertible_bond_value > 95.0);
        assert!(
            r0.cash_only_cocb_component > 0.0
                && r0.cash_only_cocb_component < r0.convertible_bond_value
        );
        assert!(r0.equity_conversion_component > 0.0);
        assert!(
            r0.effective_blended_credit_spread > 0.0
                && r0.effective_blended_credit_spread < base.credit_spread_sc
        );

        let inputs = [
            base.spot_price_s0,
            base.conversion_ratio_gamma,
            base.par_value,
            base.coupon_rate_c,
            base.risk_free_rate_r,
            base.credit_spread_sc,
            base.dividend_yield_q,
            base.equity_volatility_sigma,
            base.call_price_k,
            base.maturity_years_t,
        ];
        let mut dual_in = [Dual::<10>::from_f64(0.0); 10];
        for i in 0..10 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = TsiveriotisFernandesConvertibleSpec {
            spot_price_s0: dual_in[0],
            conversion_ratio_gamma: dual_in[1],
            par_value: dual_in[2],
            coupon_rate_c: dual_in[3],
            risk_free_rate_r: dual_in[4],
            credit_spread_sc: dual_in[5],
            dividend_yield_q: dual_in[6],
            equity_volatility_sigma: dual_in[7],
            call_price_k: dual_in[8],
            maturity_years_t: dual_in[9],
        };
        let r_ad = evaluate_tsiveriotis_fernandes_convertible_bond_credit_split(&spec_ad);
        let eps = 1e-6;
        for i in 0..10 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 10]| TsiveriotisFernandesConvertibleSpec {
                spot_price_s0: v[0],
                conversion_ratio_gamma: v[1],
                par_value: v[2],
                coupon_rate_c: v[3],
                risk_free_rate_r: v[4],
                credit_spread_sc: v[5],
                dividend_yield_q: v[6],
                equity_volatility_sigma: v[7],
                call_price_k: v[8],
                maturity_years_t: v[9],
            };
            let rp = evaluate_tsiveriotis_fernandes_convertible_bond_credit_split(&mk(p));
            let rm = evaluate_tsiveriotis_fernandes_convertible_bond_credit_split(&mk(m));
            let fd = (rp.convertible_bond_value - rm.convertible_bond_value) / (2.0 * eps);
            assert!(
                r_ad.convertible_bond_value.d[i].abs() > 1e-7,
                "zero derivative at {i}"
            );
            assert!((r_ad.convertible_bond_value.d[i] - fd).abs() < 1e-4);
        }
    }
}
