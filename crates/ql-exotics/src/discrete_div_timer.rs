use ql_ad::Real;
use ql_models::HestonParams;
use ql_vol::black_scholes_price;

#[derive(Copy, Clone, Debug)]
pub struct DiscreteDividend<T: Real = f64> {
    pub time: f64,
    pub cash: T,
    pub proportional: T,
}
#[allow(clippy::too_many_arguments)]
pub fn price_european_discrete_dividends_bv<T: Real>(
    spot: T,
    strike: T,
    vol: T,
    rate: T,
    div: T,
    expiry: T,
    dividends: &[DiscreteDividend<T>],
    is_call: bool,
) -> T {
    let mut prop = T::one();
    let mut ds = T::zero();
    let mut dk = T::zero();
    let exp_safe = expiry.max_val(T::from_f64(1e-12));
    for d in dividends
        .iter()
        .filter(|d| d.time > 0.0 && d.time <= expiry.val())
    {
        prop *= T::one() - d.proportional;
        let pv = d.cash * (-rate * d.time).exp();
        ds += pv * (expiry - T::from_f64(d.time)) / exp_safe;
        dk += pv * d.time / exp_safe;
    }
    let sstar = (spot * prop * (-div * expiry).exp() - ds) * (div * expiry).exp();
    let kstar = strike + dk * (rate * expiry).exp();
    black_scholes_price(sstar, kstar, vol, expiry, rate, div, is_call)
}
#[derive(Copy, Clone, Debug)]
pub struct TimerOptionSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub rate: T,
    pub div_yield: T,
    pub variance_budget: T,
    pub is_call: bool,
}
#[derive(Copy, Clone, Debug)]
pub struct TimerOptionValuation<T: Real> {
    pub npv: T,
    pub effective_expiry: T,
    pub variance_budget: T,
}
pub fn price_timer_option_analytic<T: Real>(
    s: &TimerOptionSpec<T>,
    h: &HestonParams<T>,
) -> TimerOptionValuation<T> {
    let b = s.variance_budget.max_val(T::from_f64(1e-14));
    let mut t = (b / h.theta.max_val(T::from_f64(1e-12))).max_val(T::from_f64(1e-8));
    for _ in 0..8 {
        let u = h.kappa * t;
        let mean = if u.val().abs() < 1e-8 {
            h.v0 * t - h.kappa * (h.v0 - h.theta) * t * t * 0.5
        } else {
            h.theta * t + (h.v0 - h.theta) * (T::one() - (-u).exp()) / h.kappa
        };
        let der = (h.v0 * (-h.kappa * t).exp() + h.theta * (T::one() - (-h.kappa * t).exp()))
            .max_val(T::from_f64(1e-10));
        t = (t - (mean - b) / der).max_val(T::from_f64(1e-8));
    }
    let bs_at = |tau: T, spot_adj: T| {
        let ts = tau.max_val(T::from_f64(1e-8));
        black_scholes_price(
            spot_adj,
            s.strike,
            b.sqrt() / ts.sqrt(),
            ts,
            s.rate,
            s.div_yield,
            s.is_call,
        )
    };
    let vend = (h.v0 * (-h.kappa * t).exp() + h.theta * (T::one() - (-h.kappa * t).exp()))
        .max_val(T::from_f64(1e-8));
    let den = h.kappa.max_val(T::from_f64(1e-8)) * b.sqrt() * vend.sqrt();
    let spot_lev = s.spot
        * (h.rho * h.sigma_v * (vend - h.v0 + h.kappa * (b - h.theta * t)) / den * 0.1).exp();
    let v0 = bs_at(t, spot_lev);
    let de = (t * 1e-3).max_val(T::from_f64(1e-5));
    let d2 = (bs_at(t + de, spot_lev) - v0 * 2.0 + bs_at(t - de, spot_lev)) / (de * de);
    let var_tau = h.sigma_v * h.sigma_v * b / (h.kappa.max_val(T::from_f64(1e-4)) * vend).powf(2.0);
    TimerOptionValuation {
        npv: (v0 + d2 * var_tau * 0.5).max_val(T::zero()),
        effective_expiry: t,
        variance_budget: b,
    }
}
#[allow(clippy::too_many_arguments)]
pub fn price_timer_option_heston_mc_ad<T: Real>(
    s: &TimerOptionSpec<T>,
    h: &HestonParams<T>,
    n_steps: usize,
    n_paths: usize,
    seed: u64,
) -> TimerOptionValuation<T> {
    let base = price_timer_option_analytic(s, h);
    let b = s.variance_budget.max_val(T::from_f64(1e-14));
    let steps = n_steps.max(8);
    let paths = n_paths.max(2);
    let horizon = base.effective_expiry * 2.5;
    let dt = horizon / steps as f64;
    let sqrt_dt = dt.sqrt();
    let rho_comp = (T::one() - h.rho * h.rho).max_val(T::zero()).sqrt();
    let rng = ql_mc::PhiloxRng::new(seed);
    let mut sum = T::zero();
    let mut tau_sum = T::zero();
    for path in 0..paths {
        let anti = path % 2 == 1;
        let id = (path / 2) as u64;
        let mut log_s = s.spot.ln();
        let mut v = h.v0.max_val(T::from_f64(1e-8));
        let mut iv = T::zero();
        let mut tau = T::zero();
        for step in 0..steps {
            let q = rng.normal4(id, step as u32);
            let sign = if anti { -1.0 } else { 1.0 };
            let z1 = q[0] * sign;
            let z2 = q[1] * sign;
            let vp = v.max_val(T::from_f64(1e-8));
            let div = vp * dt;
            if (iv + div).val() >= b.val() {
                let frac = ((b - iv) / div.max_val(T::from_f64(1e-12)))
                    .max_val(T::zero())
                    .min_val(T::one());
                let de = dt * frac;
                log_s += (s.rate - s.div_yield - vp * 0.5) * de
                    + vp.sqrt() * (h.rho * z1 + rho_comp * z2) * de.sqrt();
                tau += de;
                break;
            }
            log_s += (s.rate - s.div_yield - vp * 0.5) * dt
                + vp.sqrt() * (h.rho * z1 + rho_comp * z2) * sqrt_dt;
            v = (v + h.kappa * (h.theta - vp) * dt + h.sigma_v * vp.sqrt() * z1 * sqrt_dt)
                .max_val(T::from_f64(1e-8));
            iv += div;
            tau += dt;
        }
        let st = log_s.exp();
        let diff = if s.is_call {
            st - s.strike
        } else {
            s.strike - st
        };
        let e = 1e-4;
        sum += ((diff + (diff * diff + e * e).sqrt()) * 0.5) * (-s.rate * tau).exp();
        tau_sum += tau;
    }
    TimerOptionValuation {
        npv: sum / paths as f64,
        effective_expiry: tau_sum / paths as f64,
        variance_budget: b,
    }
}
#[derive(Copy, Clone, Debug)]
pub struct PassportOptionSpec<T: Real = f64> {
    pub spot: T,
    pub initial_account: T,
    pub rate: T,
    pub div_yield: T,
    pub account_rate: T,
    pub cost_of_carry_rate: T,
    pub vol: T,
    pub expiry: T,
    pub n_steps: usize,
}
#[allow(clippy::needless_range_loop)]
pub fn price_passport_option<T: Real>(s: &PassportOptionSpec<T>) -> T {
    let st = s.vol * s.expiry.sqrt();
    let sym = s.spot
        * (-s.div_yield * s.expiry).exp()
        * ((T::one() + st * st * 0.5) * ((st * 0.5).norm_cdf() * 2.0 - 1.0)
            + st * T::from_f64((0.5 / std::f64::consts::PI).sqrt()) * (-st * st / 8.0).exp()
            - st * st * 0.5)
        * 0.5;
    if s.initial_account.val().abs() <= 1e-12
        && (s.account_rate - s.div_yield).val().abs() <= 1e-12
        && (s.cost_of_carry_rate - s.div_yield).val().abs() <= 1e-12
    {
        return sym;
    }
    let n = s.n_steps.clamp(8, 64);
    let dt = s.expiry / n as f64;
    let step = s.vol * dt.sqrt();
    let uf = step.exp();
    let df = (-step).exp();
    let gs = ((s.rate - s.div_yield) * dt).exp();
    let prob = ((gs - df) / (uf - df).max_val(T::from_f64(1e-12)))
        .max_val(T::from_f64(0.01))
        .min_val(T::from_f64(0.99));
    let dfs = (-s.rate * dt).exp();
    let carry = ((s.rate - s.cost_of_carry_rate) * dt).exp();
    let acc = (s.account_rate * dt).exp();
    let m = 2 * n + 1;
    let wm = (s.initial_account.abs() / s.spot.max_val(T::from_f64(1e-8)) + st * 2.5)
        .max_val(T::from_f64(0.5));
    let dw = (wm * 2.0) / (m - 1) as f64;
    let grid = |i: usize| -wm + dw * i as f64;
    let smooth = |x: T| ((x + (x * x + 1e-8).sqrt()) * 0.5).max_val(T::zero());
    let interp = |arr: &[T], w: T| {
        let pos = ((w + wm) / dw)
            .max_val(T::zero())
            .min_val(T::from_f64((m - 1) as f64));
        let i = (pos.val().floor() as usize).min(m - 2);
        let a = (pos - i as f64).max_val(T::zero()).min_val(T::one());
        arr[i] * (T::one() - a) + arr[i + 1] * a
    };
    let eval_tree = |acc_step: T, carry_step: T| -> Vec<T> {
        let mut vals: Vec<T> = (0..m).map(|i| smooth(grid(i))).collect();
        for _ in 0..n {
            let mut next = vec![T::zero(); m];
            for i in 0..m {
                let w = grid(i);
                let wu_p = (w * acc_step + (uf - carry_step)) / uf;
                let wd_p = (w * acc_step + (df - carry_step)) / df;
                let cp = dfs
                    * (prob * uf * interp(&vals, wu_p)
                        + (T::one() - prob) * df * interp(&vals, wd_p));
                let wu_m = (w * acc_step - (uf - carry_step)) / uf;
                let wd_m = (w * acc_step - (df - carry_step)) / df;
                let cm = dfs
                    * (prob * uf * interp(&vals, wu_m)
                        + (T::one() - prob) * df * interp(&vals, wd_m));
                let diff = cp - cm;
                next[i] = (cp + cm + (diff * diff + 1e-12).sqrt()) * 0.5;
            }
            vals = next;
        }
        vals
    };
    let vals = eval_tree(acc, carry);
    let sym_step = (s.div_yield * dt).exp();
    let vals_sym = eval_tree(sym_step, gs);
    let w0 = s.initial_account / s.spot.max_val(T::from_f64(1e-8));
    let raw = s.spot * interp(&vals, w0);
    let raw0_sym = s.spot * interp(&vals_sym, T::zero());
    (raw - raw0_sym + sym).max_val(T::zero())
}

#[derive(Copy, Clone, Debug)]
pub struct FiniteCapTimerOptionSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub rate: T,
    pub div_yield: T,
    pub variance_budget: T,
    pub t_max: T,
    pub is_call: bool,
}

pub fn price_capped_timer_option_quadrature<T: Real>(
    s: &FiniteCapTimerOptionSpec<T>,
    h: &HestonParams<T>,
) -> TimerOptionValuation<T> {
    let base = price_timer_option_analytic(
        &TimerOptionSpec {
            spot: s.spot,
            strike: s.strike,
            rate: s.rate,
            div_yield: s.div_yield,
            variance_budget: s.variance_budget,
            is_call: s.is_call,
        },
        h,
    );
    let mu_tau = base.effective_expiry.max_val(T::from_f64(1e-6));
    let t_cap = s.t_max.max_val(T::from_f64(1e-6));
    let v_bar = (h.v0 * 0.5 + h.theta * 0.5).max_val(T::from_f64(1e-6));
    let sig_tau = (h.sigma_v / (h.kappa.max_val(T::from_f64(0.1)) * v_bar.sqrt()))
        .max_val(T::from_f64(0.05))
        * mu_tau;

    // 8-point Gauss-Legendre quadrature on [0, t_cap] for the Inverse-Gaussian stopping time density
    let gl_x = [
        -0.9602898564975363,
        -0.7966664774136267,
        -0.525_532_409_916_329,
        -0.1834346424956498,
        0.1834346424956498,
        0.525_532_409_916_329,
        0.7966664774136267,
        0.9602898564975363,
    ];
    let gl_w = [
        0.1012285362903763,
        0.2223810344533745,
        0.3137066458778873,
        0.362_683_783_378_362,
        0.362_683_783_378_362,
        0.3137066458778873,
        0.2223810344533745,
        0.1012285362903763,
    ];
    let lam_ig = (mu_tau * mu_tau * mu_tau) / (sig_tau * sig_tau).max_val(T::from_f64(1e-8));
    let b_budget = s.variance_budget.max_val(T::from_f64(1e-8));

    let mut p_hit = T::zero();
    let mut npv_hit = T::zero();
    let mut t_hit = T::zero();
    for (&xi, &wi) in gl_x.iter().zip(gl_w.iter()) {
        let t_q = t_cap * (0.5 * (xi + 1.0));
        let dt_w = t_cap * (0.5 * wi);
        let z_ig = (t_q - mu_tau) / (mu_tau * (t_q / lam_ig).sqrt()).max_val(T::from_f64(1e-8));
        let pdf_ig = (lam_ig / (t_q * t_q * t_q * (2.0 * std::f64::consts::PI)).max_val(T::from_f64(1e-12))).sqrt()
            * (-z_ig * z_ig * 0.5).exp();
        let w_q = pdf_ig * dt_w;
        let vol_eff = (b_budget / t_q.max_val(T::from_f64(1e-8))).sqrt();
        let c_q = black_scholes_price(s.spot, s.strike, vol_eff, t_q, s.rate, s.div_yield, s.is_call);
        p_hit += w_q;
        npv_hit += w_q * c_q;
        t_hit += w_q * t_q;
    }
    let p_hit_clamped = p_hit.min_val(T::one()).max_val(T::zero());
    let p_cap = T::one() - p_hit_clamped;
    let var_at_cap = (v_bar * t_cap).min_val(b_budget);
    let vol_cap = (var_at_cap / t_cap).sqrt();
    let c_cap = black_scholes_price(s.spot, s.strike, vol_cap, t_cap, s.rate, s.div_yield, s.is_call);

    let scale_hit = if p_hit.val() > 1.0 {
        T::one() / p_hit
    } else {
        T::one()
    };
    let npv = npv_hit * scale_hit + p_cap * c_cap;
    let eff_t = t_hit * scale_hit + p_cap * t_cap;
    let eff_v = b_budget * p_hit_clamped + var_at_cap * p_cap;
    TimerOptionValuation {
        npv,
        effective_expiry: eff_t,
        variance_budget: eff_v,
    }
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn discrete_parity_timer_passport() {
        let ds = vec![DiscreteDividend {
            time: 0.5,
            cash: 1.0,
            proportional: 0.01,
        }];
        let c = price_european_discrete_dividends_bv(100.0, 100.0, 0.2, 0.02, 0.01, 1.0, &ds, true);
        let p =
            price_european_discrete_dividends_bv(100.0, 100.0, 0.2, 0.02, 0.01, 1.0, &ds, false);
        let expected = 100.0 * 0.99 * (-0.01f64).exp()
            - 100.0 * (-0.02f64).exp()
            - 1.0 * (-0.02f64 * 0.5).exp();
        assert!((c - p - expected).abs() < 1e-8);
        let h_eps = 1e-5;
        let fd_exp = (price_european_discrete_dividends_bv(
            100.0,
            100.0,
            0.2,
            0.02,
            0.01,
            1.0 + h_eps,
            &ds,
            true,
        ) - price_european_discrete_dividends_bv(
            100.0,
            100.0,
            0.2,
            0.02,
            0.01,
            1.0 - h_eps,
            &ds,
            true,
        )) / (2.0 * h_eps);
        let ds_d = vec![DiscreteDividend {
            time: 0.5,
            cash: Dual::<1>::from_value(1.0),
            proportional: Dual::<1>::from_value(0.01),
        }];
        let ad_exp = price_european_discrete_dividends_bv(
            Dual::from_value(100.0),
            Dual::from_value(100.0),
            Dual::from_value(0.2),
            Dual::from_value(0.02),
            Dual::from_value(0.01),
            Dual::variable(1.0, 0),
            &ds_d,
            true,
        );
        assert!((ad_exp.d[0] - fd_exp).abs() < 1e-5);
        let h = HestonParams {
            v0: 0.04,
            kappa: 1.5,
            theta: 0.04,
            sigma_v: 0.3,
            rho: -0.5,
        };
        let t_spec = TimerOptionSpec {
            spot: 100.0,
            strike: 100.0,
            rate: 0.02,
            div_yield: 0.01,
            variance_budget: 0.04,
            is_call: true,
        };
        let t_an = price_timer_option_analytic(&t_spec, &h);
        let t_mc = price_timer_option_heston_mc_ad(&t_spec, &h, 32, 128, 42);
        assert!(t_an.npv > 0.0 && t_mc.npv > 0.0 && (t_an.npv - t_mc.npv).abs() < 2.0);
        let call_atm = black_scholes_price(100.0, 100.0, 0.2, 1.0, 0.02, 0.01, true);
        let pass0 = price_passport_option(&PassportOptionSpec {
            spot: 100.0,
            initial_account: 0.0,
            rate: 0.02,
            div_yield: 0.01,
            account_rate: 0.02,
            cost_of_carry_rate: 0.02,
            vol: 0.2,
            expiry: 1.0,
            n_steps: 16,
        });
        assert!(pass0 > 0.0 && pass0 < call_atm);
        let fd_pass = (price_passport_option(&PassportOptionSpec {
            spot: 100.0 + h_eps,
            initial_account: 2.0,
            rate: 0.02,
            div_yield: 0.01,
            account_rate: 0.02,
            cost_of_carry_rate: 0.02,
            vol: 0.2,
            expiry: 1.0,
            n_steps: 16,
        }) - price_passport_option(&PassportOptionSpec {
            spot: 100.0 - h_eps,
            initial_account: 2.0,
            rate: 0.02,
            div_yield: 0.01,
            account_rate: 0.02,
            cost_of_carry_rate: 0.02,
            vol: 0.2,
            expiry: 1.0,
            n_steps: 16,
        })) / (2.0 * h_eps);
        let d = price_passport_option(&PassportOptionSpec::<Dual<1>> {
            spot: Dual::variable(100.0, 0),
            initial_account: Dual::from_value(2.0),
            rate: Dual::from_value(0.02),
            div_yield: Dual::from_value(0.01),
            account_rate: Dual::from_value(0.02),
            cost_of_carry_rate: Dual::from_value(0.02),
            vol: Dual::from_value(0.2),
            expiry: Dual::from_value(1.0),
            n_steps: 16,
        });
        assert!((d.d[0] - fd_pass).abs() < 1e-4);
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn capped_timer_quadrature_monotonicity_limit_and_ad() {
        let h = HestonParams {
            v0: 0.04,
            kappa: 1.5,
            theta: 0.04,
            sigma_v: 0.25,
            rho: 0.0,
        };
        let s_large = FiniteCapTimerOptionSpec {
            spot: 100.0,
            strike: 100.0,
            rate: 0.0,
            div_yield: 0.0,
            variance_budget: 0.04,
            t_max: 20.0,
            is_call: true,
        };
        let mut s_tight = s_large;
        s_tight.t_max = 0.25;
        let r_large = price_capped_timer_option_quadrature(&s_large, &h);
        let r_tight = price_capped_timer_option_quadrature(&s_tight, &h);
        let r_uncapped = price_timer_option_analytic(
            &TimerOptionSpec {
                spot: 100.0,
                strike: 100.0,
                rate: 0.0,
                div_yield: 0.0,
                variance_budget: 0.04,
                is_call: true,
            },
            &h,
        );
        assert!(r_tight.npv < r_large.npv);
        assert!((r_large.npv - r_uncapped.npv).abs() < 5e-2);

        let eps = 1e-5;
        let eval_tmax = |tm: f64| {
            let mut sp = s_tight;
            sp.t_max = tm;
            price_capped_timer_option_quadrature(&sp, &h).npv
        };
        let fd = (eval_tmax(0.5 + eps) - eval_tmax(0.5 - eps)) / (2.0 * eps);
        type D = Dual<1>;
        let hd = HestonParams {
            v0: D::from_value(0.04),
            kappa: D::from_value(1.5),
            theta: D::from_value(0.04),
            sigma_v: D::from_value(0.25),
            rho: D::zero(),
        };
        let sd = FiniteCapTimerOptionSpec {
            spot: D::from_value(100.0),
            strike: D::from_value(100.0),
            rate: D::zero(),
            div_yield: D::zero(),
            variance_budget: D::from_value(0.04),
            t_max: D::variable(0.5, 0),
            is_call: true,
        };
        let ad = price_capped_timer_option_quadrature(&sd, &hd).npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
