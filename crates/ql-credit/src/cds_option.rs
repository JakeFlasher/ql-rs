#![allow(clippy::excessive_precision)]
//! Black's option on a forward CDS spread, including front-end protection.
use crate::SurvivalProbabilityCurve;
use ql_ad::Real;
use ql_curves::YieldTermStructure;

/// Black (1976) CDS option.  `risky_annuity` is the PV of one unit of
/// running spread, `default_prob_to_expiry` is the cumulative default
/// probability to option expiry, and `discount_to_expiry` is its discount
/// factor.  A payer option is a call on the spread; a receiver is a put.
#[allow(clippy::too_many_arguments)]
pub fn black_cds_option_price<T: Real>(
    forward_spread: T,
    strike_spread: T,
    risky_annuity: T,
    vol: T,
    expiry: T,
    default_prob_to_expiry: T,
    discount_to_expiry: T,
    recovery: T,
    notional: T,
    is_payer: bool,
    knocks_out: bool,
) -> T {
    let sqrt_t = expiry.max_val(T::zero()).sqrt();
    let std = vol.abs() * sqrt_t;
    let intrinsic = if is_payer {
        (forward_spread - strike_spread).max_val(T::zero())
    } else {
        (strike_spread - forward_spread).max_val(T::zero())
    };
    let option = if std.val() < 1e-14 || forward_spread.val() <= 0.0 || strike_spread.val() <= 0.0 {
        risky_annuity * intrinsic
    } else {
        let d1 = ((forward_spread / strike_spread).ln() + std * std * 0.5) / std;
        let d2 = d1 - std;
        if is_payer {
            risky_annuity * (forward_spread * d1.norm_cdf() - strike_spread * d2.norm_cdf())
        } else {
            risky_annuity * (strike_spread * (-d2).norm_cdf() - forward_spread * (-d1).norm_cdf())
        }
    } * notional;
    let fep = if !knocks_out && is_payer {
        notional * (T::one() - recovery) * default_prob_to_expiry * discount_to_expiry
    } else {
        T::zero()
    };
    option + fep
}

#[derive(Clone, Debug)]
pub struct CdsSwaptionSpec<T: Real = f64> {
    pub notional: T,
    pub strike_spread: T,
    pub option_expiry: f64,
    pub maturity: f64,
    pub recovery: T,
    pub running_spread: T,
    pub vol: T,
    pub payment_times: Vec<f64>,
    pub knocks_out_on_default: bool,
    pub is_payer: bool,
}
#[derive(Copy, Clone, Debug)]
pub struct CdsSwaptionValuation<T: Real> {
    pub npv: T,
    pub forward_spread: T,
    pub risky_annuity: T,
    pub protection_pv: T,
    pub fep: T,
}

pub fn price_cds_swaption_from_curves<T: Real>(
    s: &CdsSwaptionSpec<T>,
    discount: &impl ql_curves::YieldTermStructure<T>,
    credit: &crate::SurvivalProbabilityCurve<T>,
) -> CdsSwaptionValuation<T> {
    let te = s.option_expiry;
    let mut prot = T::zero();
    let mut ann = T::zero();
    let mut prev_surv = credit.survival_probability(te);
    let mut prev_t = te;
    for &t in &s.payment_times {
        if t <= te || t > s.maturity + 1e-12 {
            continue;
        }
        let dt = (t - prev_t).max(0.0);
        let surv = credit.survival_probability(t);
        let df = discount.discount_factor(t);
        prot += (T::one() - s.recovery) * df * (prev_surv - surv);
        ann += df * T::from_f64(dt) * (prev_surv + surv) * 0.5;
        prev_surv = surv;
        prev_t = t;
    }
    let fep = if s.knocks_out_on_default {
        T::zero()
    } else {
        (T::one() - s.recovery)
            * discount.discount_factor(te)
            * (T::one() - credit.survival_probability(te))
    };

    let ann_safe = ann.max_val(T::from_f64(1e-14));
    let eff_fwd = (prot + fep) / ann_safe;
    let npv = black_cds_option_price(
        eff_fwd,
        s.strike_spread,
        ann,
        s.vol,
        T::from_f64(te),
        T::one() - credit.survival_probability(te),
        discount.discount_factor(te),
        s.recovery,
        s.notional,
        s.is_payer,
        true,
    );
    CdsSwaptionValuation {
        npv,
        forward_spread: eff_fwd,
        risky_annuity: ann,
        protection_pv: prot,
        fep,
    }
}

#[derive(Clone, Debug)]
pub struct IndexCdsSwaptionSpec<T: Real = f64> {
    pub notional: T,
    pub index_coupon: T,
    pub strike_spread: T,
    pub option_expiry: f64,
    pub maturity: f64,
    pub payment_times: Vec<f64>,
    pub constituent_weights: Vec<f64>,
    pub constituent_recoveries: Vec<T>,
    pub vol: T,
    pub is_payer: bool,
    pub include_fep: bool,
}
#[derive(Copy, Clone, Debug)]
pub struct IndexCdsSwaptionValuation<T: Real> {
    pub npv: T,
    pub fep_adjusted_forward_spread: T,
    pub adjusted_strike_spread: T,
    pub pool_risky_annuity: T,
    pub pool_protection_pv: T,
    pub pool_fep_pv: T,
    pub expected_surviving_factor: T,
}
pub fn price_index_cds_swaption<T: Real>(
    s: &IndexCdsSwaptionSpec<T>,
    discount: &impl ql_curves::YieldTermStructure<T>,
    curves: &[crate::SurvivalProbabilityCurve<T>],
) -> IndexCdsSwaptionValuation<T> {
    let n = curves
        .len()
        .min(s.constituent_weights.len())
        .min(s.constituent_recoveries.len());
    if n == 0 {
        return IndexCdsSwaptionValuation {
            npv: T::zero(),
            fep_adjusted_forward_spread: T::zero(),
            adjusted_strike_spread: s.strike_spread,
            pool_risky_annuity: T::zero(),
            pool_protection_pv: T::zero(),
            pool_fep_pv: T::zero(),
            expected_surviving_factor: T::zero(),
        };
    }
    let wsum = s.constituent_weights[..n].iter().sum::<f64>().max(1e-14);
    let mut ann = T::zero();
    let mut prot = T::zero();
    let mut fep = T::zero();
    let mut surv = T::zero();
    let te = s.option_expiry;
    for i in 0..n {
        let w = T::from_f64(s.constituent_weights[i] / wsum);
        let rec = s.constituent_recoveries[i];
        let q_te = curves[i].survival_probability(te);
        let mut c_prot = T::zero();
        let mut c_ann = T::zero();
        let mut prev_surv = q_te;
        let mut prev_t = te;
        for &t in &s.payment_times {
            if t <= te || t > s.maturity + 1e-12 {
                continue;
            }
            let dt = (t - prev_t).max(0.0);
            let surv_t = curves[i].survival_probability(t);
            let df = discount.discount_factor(t);
            c_prot += (T::one() - rec) * df * (prev_surv - surv_t);
            c_ann += df * T::from_f64(dt) * (prev_surv + surv_t) * 0.5;
            prev_surv = surv_t;
            prev_t = t;
        }
        ann += w * c_ann;
        prot += w * c_prot;
        surv += w * q_te;
        if s.include_fep {
            fep += w * (T::one() - rec) * discount.discount_factor(te) * (T::one() - q_te);
        }
    }
    let fwd = (prot + fep) / ann.max_val(T::from_f64(1e-14));
    let kadj = s.index_coupon + (s.strike_spread - s.index_coupon) * surv;
    let npv = crate::black_cds_option_price(
        fwd,
        kadj,
        ann,
        s.vol,
        T::from_f64(s.option_expiry),
        T::one() - surv,
        discount.discount_factor(s.option_expiry),
        T::from_f64(0.4),
        s.notional,
        s.is_payer,
        true,
    );
    IndexCdsSwaptionValuation {
        npv,
        fep_adjusted_forward_spread: fwd,
        adjusted_strike_spread: kadj,
        pool_risky_annuity: ann,
        pool_protection_pv: prot,
        pool_fep_pv: fep,
        expected_surviving_factor: surv,
    }
}

#[derive(Clone, Debug)]
pub struct CmCdsMarketModelSpec<T: Real = f64> {
    pub notional: T,
    pub grid_times: Vec<f64>,
    pub start_idx: usize,
    pub end_idx: usize,
    pub cm_tenor_periods: usize,
    pub participation_rate: T,
    pub lgd: T,
    pub forward_cds_vols: Vec<T>,
    pub forward_cds_corr: T,
}

#[derive(Clone, Debug)]
pub struct CmCdsValuation<T: Real> {
    pub seller_npv: T,
    pub cmcds_premium_leg_pv: T,
    pub protection_leg_pv: T,
    pub fair_participation_rate: T,
    pub convexity_adjustment_pv: T,
    pub projected_cm_rates: Vec<T>,
}

pub fn price_cmcds_market_model<T: Real>(
    spec: &CmCdsMarketModelSpec<T>,
    discount: &impl ql_curves::YieldTermStructure<T>,
    credit: &crate::SurvivalProbabilityCurve<T>,
) -> CmCdsValuation<T> {
    let a = spec.start_idx;
    let b = spec.end_idx;
    let c = spec.cm_tenor_periods;
    let max_idx = b + c;
    if b <= a || spec.grid_times.len() <= max_idx || spec.forward_cds_vols.len() < max_idx {
        return CmCdsValuation {
            seller_npv: T::zero(),
            cmcds_premium_leg_pv: T::zero(),
            protection_leg_pv: T::zero(),
            fair_participation_rate: T::one(),
            convexity_adjustment_pv: T::zero(),
            projected_cm_rates: Vec::new(),
        };
    }

    let mut alphas = vec![0.0; max_idx + 1];
    let mut p_def = vec![T::one(); max_idx + 1];
    let mut r_one = vec![T::zero(); max_idx + 1];
    p_def[0] = discount.discount_factor(spec.grid_times[0])
        * credit.survival_probability(spec.grid_times[0]);

    for k in 1..=max_idx {
        let t_prev = spec.grid_times[k - 1];
        let t_curr = spec.grid_times[k];
        let alpha_k = (t_curr - t_prev).max(1e-8);
        alphas[k] = alpha_k;
        let q_prev = credit.survival_probability(t_prev);
        let q_curr = credit
            .survival_probability(t_curr)
            .max_val(T::from_f64(1e-12));
        p_def[k] = discount.discount_factor(t_curr) * q_curr;
        r_one[k] = spec.lgd * (q_prev / q_curr - T::one()) / alpha_k;
    }

    let mut prem_unscaled = T::zero();
    let mut prem_zero_vol = T::zero();
    let mut prot_pv = T::zero();
    let mut cm_rates = Vec::with_capacity(b - a);

    for j in (a + 1)..=b {
        let t_reset = spec.grid_times[j - 1].max(0.0);
        let mut denom = T::zero();
        for h in j..=(j + c) {
            denom += p_def[h] * alphas[h];
        }
        let denom_safe = denom.max_val(T::from_f64(1e-14));
        let mut e_cm = T::zero();
        let mut e_cm0 = T::zero();
        for i in j..=(j + c) {
            let w_ji = (p_def[i] * alphas[i]) / denom_safe;
            let sig_i = spec.forward_cds_vols[i - 1];
            let mut drift_sum = T::zero();
            for k in (j + 1)..=i {
                let sig_k = spec.forward_cds_vols[k - 1];
                let r_k = r_one[k];
                let ratio = r_k / (r_k + spec.lgd / alphas[k]).max_val(T::from_f64(1e-12));
                drift_sum += spec.forward_cds_corr * sig_k * ratio;
            }
            let conv_factor = (sig_i * drift_sum * t_reset).exp();
            e_cm += w_ji * r_one[i] * conv_factor;
            e_cm0 += w_ji * r_one[i];
        }
        cm_rates.push(e_cm);
        let annuity_j = p_def[j] * alphas[j];
        prem_unscaled += annuity_j * e_cm;
        prem_zero_vol += annuity_j * e_cm0;
        prot_pv += annuity_j * r_one[j];
    }

    let fair_phi = prot_pv / prem_unscaled.max_val(T::from_f64(1e-14));
    let prem_pv = spec.notional * spec.participation_rate * prem_unscaled;
    let prot_total = spec.notional * prot_pv;
    let conv_adj = spec.notional * (prem_unscaled - prem_zero_vol);

    CmCdsValuation {
        seller_npv: prem_pv - prot_total,
        cmcds_premium_leg_pv: prem_pv,
        protection_leg_pv: prot_total,
        fair_participation_rate: fair_phi,
        convexity_adjustment_pv: conv_adj,
        projected_cm_rates: cm_rates,
    }
}

#[derive(Clone, Debug)]
pub struct SchonbucherDefaultableRangeAccrualCapSpec<T: Real = f64> {
    pub notional: T,
    pub payment_times: Vec<T>,
    pub risk_free_rate: T,
    pub initial_hazard_rates: Vec<T>,
    pub base_recovery_rate: T,
    pub spread_volatilities: Vec<T>,
    pub spread_cap_strike: T,
    pub accrual_spread_lower: T,
    pub accrual_spread_upper: T,
    pub range_accrual_coupon_rate: T,
    pub recovery_spread_correlation: T,
}
#[derive(Clone, Debug)]
pub struct SchonbucherDefaultableRangeAccrualCapReport<T: Real = f64> {
    pub forward_credit_spreads: Vec<T>,
    pub survival_discount_factors: Vec<T>,
    pub credit_spread_cap_npv: T,
    pub credit_spread_floor_npv: T,
    pub survival_range_accrual_probabilities: Vec<T>,
    pub defaultable_range_accrual_note_npv: T,
    pub stochastic_recovery_markdown_pv: T,
}
pub fn evaluate_schonbucher_defaultable_range_accrual_and_spread_cap<T: Real>(
    s: &SchonbucherDefaultableRangeAccrualCapSpec<T>,
) -> SchonbucherDefaultableRangeAccrualCapReport<T> {
    let n = s.payment_times.len().min(s.initial_hazard_rates.len());
    let mut spreads = Vec::with_capacity(n);
    let mut dfs = Vec::with_capacity(n);
    let mut probs = Vec::with_capacity(n);
    let mut cap = T::zero();
    let mut floor = T::zero();
    let mut note = T::zero();
    let mut markdown = T::zero();
    let mut prev = T::zero();
    for i in 0..n {
        let t = s.payment_times[i];
        let dt = (t - prev).max_val(T::from_f64(1e-8));
        prev = t;
        let hz = s.initial_hazard_rates[i].max_val(T::from_f64(1e-10));
        let rec = s.base_recovery_rate.max_val(T::zero()).min_val(T::one());
        let sp = hz * (T::one() - rec).max_val(T::from_f64(1e-8));
        let df = (-s.risk_free_rate * t).exp();
        let surv = (-hz * t).exp();
        let vol = s
            .spread_volatilities
            .get(i)
            .copied()
            .unwrap_or_else(T::zero)
            .abs();
        let st = (vol * t.sqrt()).max_val(T::from_f64(1e-8));
        let d1 = ((sp / s.spread_cap_strike.max_val(T::from_f64(1e-8))).ln()
            + T::from_f64(0.5) * st * st)
            / st;
        let d2 = d1 - st;
        let ann = s.notional * df * surv * dt;
        cap += ann * (sp * d1.norm_cdf() - s.spread_cap_strike * d2.norm_cdf());
        floor += ann * (s.spread_cap_strike * (-d2).norm_cdf() - sp * (-d1).norm_cdf());
        let zhigh = ((s.accrual_spread_upper.max_val(T::from_f64(1e-8)) / sp).ln()
            + T::from_f64(0.5) * st * st)
            / st;
        let cdf_low = if s.accrual_spread_lower.val() <= 1e-10 {
            T::zero()
        } else {
            let zlow = ((s.accrual_spread_lower / sp).ln() + T::from_f64(0.5) * st * st) / st;
            zlow.norm_cdf()
        };
        let pr = (zhigh.norm_cdf() - cdf_low)
            .max_val(T::from_f64(1e-12))
            .min_val(T::one() - T::from_f64(1e-12));
        spreads.push(sp);
        dfs.push(df * surv);
        probs.push(pr);
        note += s.notional * s.range_accrual_coupon_rate * dt * df * surv * pr;
        markdown += s.notional
            * (T::one() - rec)
            * df
            * (T::one() - surv)
            * (T::one() + s.recovery_spread_correlation * vol).abs()
            * T::from_f64(0.01);
    }
    SchonbucherDefaultableRangeAccrualCapReport {
        forward_credit_spreads: spreads,
        survival_discount_factors: dfs,
        credit_spread_cap_npv: cap,
        credit_spread_floor_npv: floor,
        survival_range_accrual_probabilities: probs,
        defaultable_range_accrual_note_npv: note,
        stochastic_recovery_markdown_pv: markdown,
    }
}

#[derive(Clone, Debug)]
pub struct MoriniBrigoIndexOptionSpec<T: Real = f64> {
    pub notional: T,
    pub strike_spread: T,
    pub option_expiry: f64,
    pub payment_times: Vec<f64>,
    pub recovery_rate: T,
    pub spread_vol: T,
    pub num_constituents: usize,
    pub systemic_copula_corr: T,
    pub is_payer: bool,
}

#[derive(Copy, Clone, Debug)]
pub struct MoriniBrigoIndexOptionValuation<T: Real = f64> {
    pub npv: T,
    pub standard_market_npv: T,
    pub armageddon_probability: T,
    pub armageddon_pv_term_o2: T,
    pub loss_adjusted_forward_spread: T,
    pub no_armageddon_forward_spread: T,
    pub risky_annuity: T,
}

const GH16_X: [f64; 16] = [
    -4.688_738_939_305_818,
    -3.869_447_904_860_122_7,
    -3.176_999_161_979_956,
    -2.546_202_157_847_481_4,
    -1.951_787_990_916_254,
    -1.380_258_539_198_880_8,
    -0.822_951_449_144_655_9,
    -0.273_481_046_138_152_45,
    0.273_481_046_138_152_45,
    0.822_951_449_144_655_9,
    1.380_258_539_198_880_8,
    1.951_787_990_916_254,
    2.546_202_157_847_481_4,
    3.176_999_161_979_956,
    3.869_447_904_860_122_7,
    4.688_738_939_305_818,
];
const GH16_W: [f64; 16] = [
    2.654_807_474_011_182e-10,
    2.320_980_844_865_211e-7,
    2.711_860_092_537_881_5e-5,
    9.322_840_086_241_805e-4,
    0.012_880_311_535_509_974,
    0.083_810_041_398_985_83,
    0.280_647_458_528_533_7,
    0.507_929_479_016_613_7,
    0.507_929_479_016_613_7,
    0.280_647_458_528_533_7,
    0.083_810_041_398_985_83,
    0.012_880_311_535_509_974,
    9.322_840_086_241_805e-4,
    2.711_860_092_537_881_5e-5,
    2.320_980_844_865_211e-7,
    2.654_807_474_011_182e-10,
];

fn inv_norm_cdf_generic<T: Real>(p: T) -> T {
    let p_clamped = p.val().clamp(1e-12, 1.0 - 1e-12);
    let x0 = ql_math::inv_norm_cdf(p_clamped);
    let mut x = T::from_f64(x0);
    for _ in 0..3 {
        let err = x.norm_cdf() - p;
        let pdf = x.norm_pdf().max_val(T::from_f64(1e-14));
        x -= err / pdf;
    }
    x
}

/// Morini & Brigo (2008 / 2011) Arbitrage-Free "No-Armageddon" Credit Index Option engine
/// (Eqs. 3-4, 8, 21, 22, 25).
pub fn price_morini_brigo_credit_index_option<T: Real>(
    spec: &MoriniBrigoIndexOptionSpec<T>,
    discount_curve: &impl YieldTermStructure<T>,
    survival_curve: &SurvivalProbabilityCurve<T>,
) -> MoriniBrigoIndexOptionValuation<T> {
    let ta = spec.option_expiry.max(1e-6);
    let lgd = (T::one() - spec.recovery_rate).max_val(T::from_f64(1e-8));
    let n_names = spec.num_constituents.max(1);
    let rho = spec
        .systemic_copula_corr
        .max_val(T::zero())
        .min_val(T::from_f64(0.999_9));

    let mut annuity = T::zero();
    let mut prot_leg = T::zero();
    let mut prev_t = ta;
    let mut prev_q = survival_curve.survival_probability(ta);
    for &tj in &spec.payment_times {
        if tj <= ta + 1e-8 {
            continue;
        }
        let dt = tj - prev_t;
        let df = discount_curve.discount_factor(tj);
        let qj = survival_curve.survival_probability(tj);
        annuity += df * qj * dt;
        prot_leg += df * lgd * (prev_q - qj);
        prev_t = tj;
        prev_q = qj;
    }
    let annuity_safe = annuity.max_val(T::from_f64(1e-12));
    let df_ta = discount_curve.discount_factor(ta);
    let q_ta = survival_curve.survival_probability(ta);
    let def_prob_ta = (T::one() - q_ta)
        .max_val(T::from_f64(1e-12))
        .min_val(T::from_f64(1.0 - 1e-12));
    let fep_pv = lgd * df_ta * def_prob_ta;
    let s_tilde = (prot_leg + fep_pv) / annuity_safe;

    // 1-Factor Gaussian copula Armageddon probability P(tau_hat <= T_A) via 16-pt Gauss-Hermite
    let c_thresh = inv_norm_cdf_generic(def_prob_ta);
    let sqrt_rho = (rho + T::from_f64(1e-16)).sqrt();
    let sqrt_1m_rho = (T::one() - rho).max_val(T::from_f64(1e-12)).sqrt();
    let inv_sqrt_pi = 1.0 / std::f64::consts::PI.sqrt();
    let mut p_arm = T::zero();
    for m in 0..16 {
        let z = T::from_f64(GH16_X[m] * std::f64::consts::SQRT_2);
        let cond_p = ((c_thresh - sqrt_rho * z) / sqrt_1m_rho)
            .norm_cdf()
            .max_val(T::from_f64(1e-15));
        let prod_n = (cond_p.ln() * (n_names as f64)).exp();
        p_arm += prod_n * (GH16_W[m] * inv_sqrt_pi);
    }
    p_arm = p_arm.min_val(def_prob_ta).max_val(T::zero());

    let o2_arm = lgd * df_ta * p_arm;
    let s_hat = ((prot_leg + fep_pv - o2_arm) / annuity_safe).max_val(T::from_f64(1e-10));
    let k = spec.strike_spread.max_val(T::from_f64(1e-10));
    let st = (spec.spread_vol * ta.sqrt()).max_val(T::from_f64(1e-8));

    let d1_hat = (s_hat / k).ln() / st + st * 0.5;
    let d2_hat = d1_hat - st;
    let npv = if spec.is_payer {
        spec.notional
            * (annuity_safe * (s_hat * d1_hat.norm_cdf() - k * d2_hat.norm_cdf()) + o2_arm)
    } else {
        spec.notional * annuity_safe * (k * (-d2_hat).norm_cdf() - s_hat * (-d1_hat).norm_cdf())
    };

    let d1_std = (s_tilde / k).ln() / st + st * 0.5;
    let d2_std = d1_std - st;
    let standard_market_npv = if spec.is_payer {
        spec.notional * annuity_safe * (s_tilde * d1_std.norm_cdf() - k * d2_std.norm_cdf())
    } else {
        spec.notional * annuity_safe * (k * (-d2_std).norm_cdf() - s_tilde * (-d1_std).norm_cdf())
    };

    MoriniBrigoIndexOptionValuation {
        npv,
        standard_market_npv,
        armageddon_probability: p_arm,
        armageddon_pv_term_o2: o2_arm,
        loss_adjusted_forward_spread: s_tilde,
        no_armageddon_forward_spread: s_hat,
        risky_annuity: annuity_safe,
    }
}

#[cfg(test)]
mod domain_suite_26_tests {
    use super::*;
    use crate::SurvivalProbabilityCurve;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn suite26_cds_swaption_curve_backed() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let disc = FlatForwardCurve::new(d, 0.02);
        let credit = SurvivalProbabilityCurve::new(d, vec![1.0, 2.0, 5.0], vec![0.02, 0.02, 0.02]);
        let mut s = CdsSwaptionSpec {
            notional: 100.0,
            strike_spread: 0.012,
            option_expiry: 1.0,
            maturity: 5.0,
            recovery: 0.4,
            running_spread: 0.012,
            vol: 0.25,
            payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            knocks_out_on_default: false,
            is_payer: true,
        };
        let pay = price_cds_swaption_from_curves(&s, &disc, &credit);
        s.is_payer = false;
        let rec = price_cds_swaption_from_curves(&s, &disc, &credit);
        let expected_diff =
            s.notional * (pay.protection_pv + pay.fep - s.strike_spread * pay.risky_annuity);
        assert!((pay.npv - rec.npv - expected_diff).abs() < 1e-10);

        let eval_vol = |v: f64| {
            let mut spec = s.clone();
            spec.is_payer = true;
            spec.vol = v;
            price_cds_swaption_from_curves(&spec, &disc, &credit).npv
        };
        let h = 1e-5;
        let fd_vega = (eval_vol(0.25 + h) - eval_vol(0.25 - h)) / (2.0 * h);
        type D = Dual<1>;
        let disc_d = FlatForwardCurve::new(d, D::from_value(0.02));
        let credit_d = SurvivalProbabilityCurve::new(
            d,
            vec![1.0, 2.0, 5.0],
            vec![
                D::from_value(0.02),
                D::from_value(0.02),
                D::from_value(0.02),
            ],
        );
        let sd = CdsSwaptionSpec {
            notional: D::from_value(100.0),
            strike_spread: D::from_value(0.012),
            option_expiry: 1.0,
            maturity: 5.0,
            recovery: D::from_value(0.4),
            running_spread: D::from_value(0.012),
            vol: D::variable(0.25, 0),
            payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            knocks_out_on_default: false,
            is_payer: true,
        };
        let ad = price_cds_swaption_from_curves(&sd, &disc_d, &credit_d).npv;
        assert!((ad.d[0] - fd_vega).abs() < 1e-6);
    }
}

#[cfg(test)]
mod domain_suite_27_tests {
    use super::*;
    use crate::SurvivalProbabilityCurve;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn index_cds_swaption_pool_parity() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let disc = FlatForwardCurve::new(d, 0.02);
        let curve = SurvivalProbabilityCurve::new(d, vec![1.0, 2.0, 5.0], vec![0.03, 0.03, 0.03]);
        let mut s = IndexCdsSwaptionSpec {
            notional: 100.0,
            index_coupon: 0.012,
            strike_spread: 0.012,
            option_expiry: 1.0,
            maturity: 5.0,
            payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            constituent_weights: vec![1.0],
            constituent_recoveries: vec![0.4],
            vol: 0.25,
            is_payer: true,
            include_fep: true,
        };
        let p = price_index_cds_swaption(&s, &disc, std::slice::from_ref(&curve));
        let single = price_cds_swaption_from_curves(
            &CdsSwaptionSpec {
                notional: 100.0,
                strike_spread: 0.012,
                option_expiry: 1.0,
                maturity: 5.0,
                recovery: 0.4,
                running_spread: 0.012,
                vol: 0.25,
                payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
                knocks_out_on_default: false,
                is_payer: true,
            },
            &disc,
            &curve,
        );
        assert!((p.npv - single.npv).abs() < 1e-10);
        s.is_payer = false;
        let r = price_index_cds_swaption(&s, &disc, std::slice::from_ref(&curve));
        assert!(
            (p.npv
                - r.npv
                - s.notional
                    * p.pool_risky_annuity
                    * (p.fep_adjusted_forward_spread - p.adjusted_strike_spread))
                .abs()
                < 1e-8
        );
        let h = 1e-5;
        let eval_rec = |rec_val: f64| {
            let mut spec = s.clone();
            spec.is_payer = true;
            spec.constituent_recoveries = vec![rec_val];
            price_index_cds_swaption(&spec, &disc, std::slice::from_ref(&curve)).npv
        };
        let fd_rec = (eval_rec(0.4 + h) - eval_rec(0.4 - h)) / (2.0 * h);
        type D = Dual<1>;
        let disc_d = FlatForwardCurve::new(d, D::from_value(0.02));
        let curve_d =
            SurvivalProbabilityCurve::new(d, vec![1.0, 2.0, 5.0], vec![D::from_value(0.03); 3]);
        let sd = IndexCdsSwaptionSpec {
            notional: D::from_value(100.0),
            index_coupon: D::from_value(0.012),
            strike_spread: D::from_value(0.012),
            option_expiry: 1.0,
            maturity: 5.0,
            payment_times: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            constituent_weights: vec![1.0],
            constituent_recoveries: vec![D::variable(0.4, 0)],
            vol: D::from_value(0.25),
            is_payer: true,
            include_fep: true,
        };
        let ad = price_index_cds_swaption(&sd, &disc_d, std::slice::from_ref(&curve_d)).npv;
        assert!((ad.d[0] - fd_rec).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_31_tests {
    use super::*;
    use crate::SurvivalProbabilityCurve;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;

    #[test]
    fn brigo_cmcds_one_period_collapse_convexity_and_ad() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let disc = FlatForwardCurve::new(d, 0.03);
        let credit = SurvivalProbabilityCurve::new(
            d,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            vec![0.02; 8],
        );
        let grid = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        // Brigo (2008) Prop 2.4 invariant: c = 0 (1-period reset) => zero convexity, fair_participation_rate == 1.0, seller_npv == 0.0
        let spec_c0 = CmCdsMarketModelSpec {
            notional: 100.0,
            grid_times: grid.clone(),
            start_idx: 0,
            end_idx: 4,
            cm_tenor_periods: 0,
            participation_rate: 1.0,
            lgd: 0.6,
            forward_cds_vols: vec![0.35; 8],
            forward_cds_corr: 0.85,
        };
        let v0 = price_cmcds_market_model(&spec_c0, &disc, &credit);
        assert!(v0.seller_npv.abs() < 1e-12);
        assert!((v0.fair_participation_rate - 1.0).abs() < 1e-12);
        assert!(v0.convexity_adjustment_pv.abs() < 1e-12);

        // c = 3 (forward 4-period CMCDS): positive convexity adjustment & fair participation rate < 1
        let mut spec_c3 = spec_c0.clone();
        spec_c3.cm_tenor_periods = 3;
        let v3 = price_cmcds_market_model(&spec_c3, &disc, &credit);
        assert!(v3.convexity_adjustment_pv > 0.0);
        assert!(v3.fair_participation_rate < 1.0);

        let h = 1e-5;
        let eval_vol = |v: f64| {
            let mut s = spec_c3.clone();
            s.forward_cds_vols = vec![v; 8];
            price_cmcds_market_model(&s, &disc, &credit).fair_participation_rate
        };
        let fd = (eval_vol(0.35 + h) - eval_vol(0.35 - h)) / (2.0 * h);
        type D = Dual<1>;
        let disc_d = FlatForwardCurve::new(d, D::from_value(0.03));
        let credit_d = SurvivalProbabilityCurve::new(
            d,
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            vec![D::from_value(0.02); 8],
        );
        let spec_d = CmCdsMarketModelSpec {
            notional: D::from_value(100.0),
            grid_times: grid,
            start_idx: 0,
            end_idx: 4,
            cm_tenor_periods: 3,
            participation_rate: D::one(),
            lgd: D::from_value(0.6),
            forward_cds_vols: vec![D::variable(0.35, 0); 8],
            forward_cds_corr: D::from_value(0.85),
        };
        let ad = price_cmcds_market_model(&spec_d, &disc_d, &credit_d).fair_participation_rate;
        assert!(ad.d[0] < 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

#[cfg(test)]
mod domain_suite_32_tests {
    use super::*;
    use crate::SurvivalProbabilityCurve;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;

    #[test]
    fn morini_brigo_no_armageddon_put_call_parity_copula_limits_and_ad() {
        let d = Date::from_ymd(2024, 1, 1).unwrap();
        let disc = FlatForwardCurve::new(d, 0.03);
        let surv = SurvivalProbabilityCurve::new(d, vec![1.0, 2.0, 3.0, 4.0, 5.0], vec![0.025; 5]);
        let payer_spec = MoriniBrigoIndexOptionSpec {
            notional: 100.0,
            strike_spread: 0.025,
            option_expiry: 1.0,
            payment_times: vec![2.0, 3.0, 4.0, 5.0],
            recovery_rate: 0.40,
            spread_vol: 0.45,
            num_constituents: 50,
            systemic_copula_corr: 0.65,
            is_payer: true,
        };
        let rec_spec = MoriniBrigoIndexOptionSpec {
            is_payer: false,
            ..payer_spec.clone()
        };
        let payer = price_morini_brigo_credit_index_option(&payer_spec, &disc, &surv);
        let rec = price_morini_brigo_credit_index_option(&rec_spec, &disc, &surv);
        assert!(payer.armageddon_probability > 0.0);
        assert!(payer.no_armageddon_forward_spread < payer.loss_adjusted_forward_spread);
        // Exact Morini-Brigo (2008) Eq. 25 Put-Call Parity:
        // Payer - Receiver == Notional * Annuity * (S_tilde - K)
        let fwd_pv = payer_spec.notional
            * payer.risky_annuity
            * (payer.loss_adjusted_forward_spread - payer_spec.strike_spread);
        assert!((payer.npv - rec.npv - fwd_pv).abs() < 1e-10);

        let h = 1e-5;
        let eval_rho = |rho: f64| {
            let mut s = payer_spec.clone();
            s.systemic_copula_corr = rho;
            price_morini_brigo_credit_index_option(&s, &disc, &surv).npv
        };
        let fd = (eval_rho(0.65 + h) - eval_rho(0.65 - h)) / (2.0 * h);
        type D = Dual<1>;
        let disc_d = FlatForwardCurve::new(d, D::from_value(0.03));
        let surv_d = SurvivalProbabilityCurve::new(
            d,
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
            vec![D::from_value(0.025); 5],
        );
        let spec_d = MoriniBrigoIndexOptionSpec {
            notional: D::from_value(100.0),
            strike_spread: D::from_value(0.025),
            option_expiry: 1.0,
            payment_times: vec![2.0, 3.0, 4.0, 5.0],
            recovery_rate: D::from_value(0.40),
            spread_vol: D::from_value(0.45),
            num_constituents: 50,
            systemic_copula_corr: D::variable(0.65, 0),
            is_payer: true,
        };
        let ad = price_morini_brigo_credit_index_option(&spec_d, &disc_d, &surv_d).npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_38_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(notional: T) -> SchonbucherDefaultableRangeAccrualCapSpec<T> {
        SchonbucherDefaultableRangeAccrualCapSpec {
            notional,
            payment_times: vec![T::from_f64(0.5), T::one(), T::from_f64(1.5)],
            risk_free_rate: T::from_f64(0.02),
            initial_hazard_rates: vec![T::from_f64(0.02); 3],
            base_recovery_rate: T::from_f64(0.4),
            spread_volatilities: vec![T::from_f64(0.3); 3],
            spread_cap_strike: T::from_f64(0.03),
            accrual_spread_lower: T::from_f64(0.0),
            accrual_spread_upper: T::from_f64(0.06),
            range_accrual_coupon_rate: T::from_f64(0.02),
            recovery_spread_correlation: T::from_f64(-0.2),
        }
    }
    #[test]
    fn schonbucher_invariants_and_dual_fd() {
        let r = evaluate_schonbucher_defaultable_range_accrual_and_spread_cap(&spec(100.0));
        assert!(r.credit_spread_cap_npv >= 0.0 && r.credit_spread_floor_npv >= 0.0);
        assert!(
            r.survival_range_accrual_probabilities
                .iter()
                .all(|p| *p > 0.0 && *p < 1.0)
        );
        let h = 1e-5;
        let fd = (evaluate_schonbucher_defaultable_range_accrual_and_spread_cap(&spec(100.0 + h))
            .defaultable_range_accrual_note_npv
            - evaluate_schonbucher_defaultable_range_accrual_and_spread_cap(&spec(100.0 - h))
                .defaultable_range_accrual_note_npv)
            / (2.0 * h);
        let ad = evaluate_schonbucher_defaultable_range_accrual_and_spread_cap(&spec(
            Dual::<1>::variable(100.0, 0),
        ))
        .defaultable_range_accrual_note_npv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

mod evaluate_brigo_morini_cmcds_and_credit_index_option_module {
    use ql_ad::Real;
    #[derive(Clone, Debug)]
    pub struct BrigoMoriniCmcdsCreditOptionSpec<T: Real = f64> {
        pub reset_times: Vec<T>,
        pub forward_cds_spreads: Vec<T>,
        pub spread_volatility: T,
        pub annuity_convexity_gamma: T,
        pub flat_par_cds_spread: T,
        pub option_expiry: T,
        pub option_strike_spread: T,
        pub risky_pv01_annuity: T,
        pub hazard_rate: T,
        pub recovery_rate: T,
        pub risk_free_rate: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct BrigoMoriniCmcdsCreditOptionReport<T: Real = f64> {
        pub cmcds_convexity_uplift_bps: T,
        pub cmcds_equilibrium_participation_rate: T,
        pub front_end_protection_pv: T,
        pub fep_adjusted_forward_spread: T,
        pub fep_adjusted_strike_spread: T,
        pub morini_brigo_payer_swaption_pv: T,
        pub naive_unadjusted_swaption_pv: T,
    }
    pub fn evaluate_brigo_morini_cmcds_and_credit_index_option<T: Real>(
        s: &BrigoMoriniCmcdsCreditOptionSpec<T>,
    ) -> BrigoMoriniCmcdsCreditOptionReport<T> {
        let mut den = T::zero();
        let mut num = T::zero();
        let mut up = T::zero();
        for i in 0..s.reset_times.len().min(s.forward_cds_spreads.len()) {
            let t = s.reset_times[i];
            let dt = if i == 0 { t } else { t - s.reset_times[i - 1] };
            let w = dt * (-s.risk_free_rate * t).exp() * (-s.hazard_rate * t).exp();
            let u =
                (s.annuity_convexity_gamma * s.spread_volatility * s.spread_volatility * t).exp();
            den += w * s.forward_cds_spreads[i] * u;
            num += w * s.flat_par_cds_spread;
            up += w * s.forward_cds_spreads[i] * (u - T::one());
        }
        let part = (num / den.max_val(T::from_f64(1e-8))).min_val(T::one());
        let df = (-s.risk_free_rate * s.option_expiry).exp();
        let surv = (-s.hazard_rate * s.option_expiry).exp();
        let fep = (T::one() - s.recovery_rate) * (T::one() - surv) * df;
        let adjf = s.flat_par_cds_spread + fep / s.risky_pv01_annuity.max_val(T::from_f64(1e-8));
        let adjk = s.option_strike_spread + fep / s.risky_pv01_annuity.max_val(T::from_f64(1e-8));
        let vol = s.spread_volatility * s.option_expiry.sqrt();
        let d1 = ((adjf / adjk.max_val(T::from_f64(1e-8))).ln() + T::from_f64(0.5) * vol * vol)
            / vol.max_val(T::from_f64(1e-8));
        let d2 = d1 - vol;
        let pv = s.risky_pv01_annuity * (adjf * d1.norm_cdf() - adjk * d2.norm_cdf()) * df;
        let d1_naive =
            ((s.flat_par_cds_spread / s.option_strike_spread.max_val(T::from_f64(1e-8))).ln()
                + T::from_f64(0.5) * vol * vol)
                / vol.max_val(T::from_f64(1e-8));
        let d2_naive = d1_naive - vol;
        let naive = s.risky_pv01_annuity
            * (s.flat_par_cds_spread * d1_naive.norm_cdf()
                - s.option_strike_spread * d2_naive.norm_cdf())
            * df;
        BrigoMoriniCmcdsCreditOptionReport {
            cmcds_convexity_uplift_bps: up / den.max_val(T::from_f64(1e-8)) * T::from_f64(1e4),
            cmcds_equilibrium_participation_rate: part,
            front_end_protection_pv: fep,
            fep_adjusted_forward_spread: adjf,
            fep_adjusted_strike_spread: adjk,
            morini_brigo_payer_swaption_pv: pv,
            naive_unadjusted_swaption_pv: naive,
        }
    }
    #[cfg(test)]
    mod domain_suite_44_tests {
        use super::*;
        use ql_ad::Dual;

        fn spec<T: Real>(x: T) -> BrigoMoriniCmcdsCreditOptionSpec<T> {
            BrigoMoriniCmcdsCreditOptionSpec {
                reset_times: vec![T::one(), T::from_f64(2.0)],
                forward_cds_spreads: vec![T::from_f64(0.01), T::from_f64(0.01)],
                spread_volatility: T::from_f64(0.2),
                annuity_convexity_gamma: T::from_f64(0.5),
                flat_par_cds_spread: T::from_f64(0.01) + x,
                option_expiry: T::one(),
                option_strike_spread: T::from_f64(0.01),
                risky_pv01_annuity: T::from_f64(4.0),
                hazard_rate: T::from_f64(0.02),
                recovery_rate: T::from_f64(0.4),
                risk_free_rate: T::from_f64(0.02),
            }
        }

        #[test]
        fn fep_positive() {
            let r = evaluate_brigo_morini_cmcds_and_credit_index_option(&spec(0.0));
            assert!(r.cmcds_equilibrium_participation_rate < 1.0);
            assert!(r.fep_adjusted_forward_spread > 0.01);
            let h = 1e-6;
            let fd = (evaluate_brigo_morini_cmcds_and_credit_index_option(&spec(h))
                .morini_brigo_payer_swaption_pv
                - evaluate_brigo_morini_cmcds_and_credit_index_option(&spec(-h))
                    .morini_brigo_payer_swaption_pv)
                / (2.0 * h);
            let ad = evaluate_brigo_morini_cmcds_and_credit_index_option(&spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .morini_brigo_payer_swaption_pv;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_brigo_morini_cmcds_and_credit_index_option_module::*;

mod evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva_module {
    use ql_ad::Real;
    #[derive(Clone, Debug)]
    pub struct CrepeyBurgardKjaerBilateralFvaSpec<T: Real = f64> {
        pub time_grid: Vec<T>,
        pub epe_profile: Vec<T>,
        pub ene_profile: Vec<T>,
        pub cpty_hazard_rate: T,
        pub bank_hazard_rate: T,
        pub cpty_recovery_rate: T,
        pub bank_recovery_rate: T,
        pub risk_free_rate: T,
        pub liquidity_funding_spread: T,
        pub burgard_kjaer_alpha: T,
    }
    #[derive(Clone, Copy, Debug)]
    pub struct CrepeyBurgardKjaerBilateralFvaReport<T: Real = f64> {
        pub bilateral_cva: T,
        pub bilateral_dva: T,
        pub burgard_kjaer_effective_funding_spread: T,
        pub funding_cost_adjustment_fca: T,
        pub funding_benefit_adjustment_fba: T,
        pub dva_fca_double_counting_overlap: T,
        pub total_burgard_kjaer_xva_charge: T,
    }
    pub fn evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva<T: Real>(
        s: &CrepeyBurgardKjaerBilateralFvaSpec<T>,
    ) -> CrepeyBurgardKjaerBilateralFvaReport<T> {
        let mut c = T::zero();
        let mut d = T::zero();
        let mut f = T::zero();
        let mut b = T::zero();
        for i in 0..s
            .time_grid
            .len()
            .min(s.epe_profile.len())
            .min(s.ene_profile.len())
        {
            let t = s.time_grid[i];
            let dt = if i == 0 { t } else { t - s.time_grid[i - 1] };
            let w = dt * (-(s.risk_free_rate + s.bank_hazard_rate + s.cpty_hazard_rate) * t).exp();
            c += w * (T::one() - s.cpty_recovery_rate) * s.cpty_hazard_rate * s.epe_profile[i];
            d += w * (T::one() - s.bank_recovery_rate) * s.bank_hazard_rate * s.ene_profile[i];
            f += w
                * (s.liquidity_funding_spread
                    + s.burgard_kjaer_alpha
                        * (T::one() - s.bank_recovery_rate)
                        * s.bank_hazard_rate)
                * s.epe_profile[i];
            b += w * s.liquidity_funding_spread * s.ene_profile[i];
        }
        let eff = s.liquidity_funding_spread
            + s.burgard_kjaer_alpha * (T::one() - s.bank_recovery_rate) * s.bank_hazard_rate;
        CrepeyBurgardKjaerBilateralFvaReport {
            bilateral_cva: c,
            bilateral_dva: d,
            burgard_kjaer_effective_funding_spread: eff,
            funding_cost_adjustment_fca: f,
            funding_benefit_adjustment_fba: b,
            dva_fca_double_counting_overlap: s.burgard_kjaer_alpha * d,
            total_burgard_kjaer_xva_charge: c - d + f - b,
        }
    }

    #[cfg(test)]
    mod domain_suite_46_tests {
        use super::*;
        use ql_ad::Dual;
        fn spec<T: Real>(x: T) -> CrepeyBurgardKjaerBilateralFvaSpec<T> {
            CrepeyBurgardKjaerBilateralFvaSpec {
                time_grid: vec![T::from_f64(0.5), T::one(), T::from_f64(2.0)],
                epe_profile: vec![T::from_f64(10.0) + x, T::from_f64(12.0), T::from_f64(8.0)],
                ene_profile: vec![T::from_f64(6.0), T::from_f64(7.0), T::from_f64(5.0)],
                cpty_hazard_rate: T::from_f64(0.02),
                bank_hazard_rate: T::from_f64(0.015),
                cpty_recovery_rate: T::from_f64(0.4),
                bank_recovery_rate: T::from_f64(0.4),
                risk_free_rate: T::from_f64(0.03),
                liquidity_funding_spread: T::from_f64(0.008),
                burgard_kjaer_alpha: T::from_f64(0.5),
            }
        }
        #[test]
        fn suite46_credit_invariants_and_dual_fd() {
            let r = evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva(&spec(0.0));
            assert!(r.bilateral_cva > 0.0 && r.bilateral_dva > 0.0);
            let h = 1e-6;
            let fd = (evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva(&spec(h))
                .total_burgard_kjaer_xva_charge
                - evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva(&spec(-h))
                    .total_burgard_kjaer_xva_charge)
                / (2.0 * h);
            let ad = evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva(&spec(
                Dual::<1>::variable(0.0, 0),
            ))
            .total_burgard_kjaer_xva_charge;
            assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_crepey_burgard_kjaer_bilateral_xva_and_fva_module::*;

mod evaluate_schonbucher_cds_swaption_and_fep_module {
    // Schönbucher (2003) & Lando (2004) Credit-Default Swaption (`Payer`/`Receiver`)
    // Risky-Annuity Numéraire & Front-End Protection (`FEP`) Black-76 Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct SchonbucherCdsSwaptionFepSpec<T: Real = f64> {
        pub hazard_rate_lambda: T,
        pub recovery_rate: T,
        pub risk_free_rate: T,
        pub option_expiry_te: T,
        pub swap_tenor_tm: T,
        pub strike_spread_k: T,
        pub spread_lognormal_vol_sigma: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct SchonbucherCdsSwaptionFepReport<T: Real = f64> {
        pub forward_risky_annuity_rpv01: T,
        pub knockout_forward_cds_spread: T,
        pub front_end_protection_pv: T,
        pub fep_adjusted_forward_spread: T,
        pub knockout_payer_swaption_bps: T,
        pub non_knockout_payer_swaption_bps: T,
        pub receiver_swaption_bps: T,
        pub risky_put_call_parity_error: T,
    }

    pub fn evaluate_schonbucher_cds_swaption_and_fep<T: Real>(
        spec: &SchonbucherCdsSwaptionFepSpec<T>,
    ) -> SchonbucherCdsSwaptionFepReport<T> {
        let lam = spec.hazard_rate_lambda.max_val(T::from_f64(1e-5));
        let r = spec.risk_free_rate.max_val(T::from_f64(1e-4));
        let lgd = T::one() - spec.recovery_rate;
        let te = spec.option_expiry_te.max_val(T::from_f64(0.05));
        let tm = spec.swap_tenor_tm.max_val(T::from_f64(0.5));
        let k = spec.strike_spread_k.max_val(T::from_f64(1e-5));
        let sig = spec.spread_lognormal_vol_sigma.max_val(T::from_f64(1e-4));

        // Forward survival-contingent risky annuity A(0; T_e, T_e + T_M):
        let u = r + lam;
        let df_surv_te = (-u * te).exp();
        let fwd_rpv01 = df_surv_te * (T::one() - (-u * tm).exp()) / u;

        // Knockout forward CDS par spread S_fwd = (1 - R) * lambda:
        let s_fwd = lgd * lam;

        // Front-End Protection (FEP) PV for default prior to option expiry T_e:
        let fep_pv = lgd * lam * (T::one() - df_surv_te) / u;
        let s_fep = s_fwd + fep_pv / fwd_rpv01.max_val(T::from_f64(1e-8));

        // Black-76 under risky annuity numéraire for knockout payer/receiver swaption:
        let std = sig * te.sqrt();
        let d1_ko = ((s_fwd / k).ln() + T::from_f64(0.5) * std * std) / std;
        let d2_ko = d1_ko - std;
        let payer_ko = fwd_rpv01 * (s_fwd * norm_cdf(d1_ko) - k * norm_cdf(d2_ko));
        let receiver_ko = fwd_rpv01 * (k * norm_cdf(-d2_ko) - s_fwd * norm_cdf(-d1_ko));

        // Non-knockout payer swaption using FEP-adjusted forward spread S_FEP:
        let d1_fep = ((s_fep / k).ln() + T::from_f64(0.5) * std * std) / std;
        let d2_fep = d1_fep - std;
        let payer_fep = fwd_rpv01 * (s_fep * norm_cdf(d1_fep) - k * norm_cdf(d2_fep));

        // Risky-annuity Put-Call Parity: Payer_KO - Receiver_KO == A(0; T_e, T_M) * (S_fwd - K)
        let pcp_err = ((payer_ko - receiver_ko) - fwd_rpv01 * (s_fwd - k)).abs();

        SchonbucherCdsSwaptionFepReport {
            forward_risky_annuity_rpv01: fwd_rpv01,
            knockout_forward_cds_spread: s_fwd,
            front_end_protection_pv: fep_pv,
            fep_adjusted_forward_spread: s_fep,
            knockout_payer_swaption_bps: payer_ko * T::from_f64(1e4),
            non_knockout_payer_swaption_bps: payer_fep * T::from_f64(1e4),
            receiver_swaption_bps: receiver_ko * T::from_f64(1e4),
            risky_put_call_parity_error: pcp_err,
        }
    }

    #[cfg(test)]
    mod domain_suite_52_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn schonbucher_cds_swaption_fep_invariants_and_ad_parity() {
            let spec = SchonbucherCdsSwaptionFepSpec {
                hazard_rate_lambda: 0.022,
                recovery_rate: 0.40,
                risk_free_rate: 0.035,
                option_expiry_te: 1.0,
                swap_tenor_tm: 5.0,
                strike_spread_k: 0.013,
                spread_lognormal_vol_sigma: 0.45,
            };
            let rep = evaluate_schonbucher_cds_swaption_and_fep(&spec);
            assert!(rep.risky_put_call_parity_error < 1e-12);
            assert!(rep.non_knockout_payer_swaption_bps > rep.knockout_payer_swaption_bps);

            let h = 1e-7;
            let r_ad = evaluate_schonbucher_cds_swaption_and_fep(&SchonbucherCdsSwaptionFepSpec {
                hazard_rate_lambda: Dual::<1>::variable(spec.hazard_rate_lambda, 0),
                recovery_rate: Dual::<1>::from_value(spec.recovery_rate),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                option_expiry_te: Dual::<1>::from_value(spec.option_expiry_te),
                swap_tenor_tm: Dual::<1>::from_value(spec.swap_tenor_tm),
                strike_spread_k: Dual::<1>::from_value(spec.strike_spread_k),
                spread_lognormal_vol_sigma: Dual::<1>::from_value(spec.spread_lognormal_vol_sigma),
            });
            let mut sp = spec;
            sp.hazard_rate_lambda += h;
            let mut sm = spec;
            sm.hazard_rate_lambda -= h;
            let fd = (evaluate_schonbucher_cds_swaption_and_fep(&sp)
                .non_knockout_payer_swaption_bps
                - evaluate_schonbucher_cds_swaption_and_fep(&sm).non_knockout_payer_swaption_bps)
                / (2.0 * h);
            assert!((r_ad.non_knockout_payer_swaption_bps.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_schonbucher_cds_swaption_and_fep_module::*;

mod evaluate_okane_isda_standard_cds_upfront_triangle_module {
    // O'Kane (2008) & Brigo–Mercurio (2006) ISDA Standard CDS Upfront,
    // Dirty RPV01 & CS01/IR01/JTD Triangle Hedge Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct OKaneIsdaStandardCdsUpfrontSpec<T: Real = f64> {
        pub notional_amount: T,
        pub market_par_spread_s: T,
        pub standard_running_coupon_c: T,
        pub recovery_rate_r: T,
        pub ois_zero_rate_r: T,
        pub maturity_years_t: T,
        pub accrued_fraction_alpha: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct OKaneIsdaStandardCdsUpfrontReport<T: Real = f64> {
        pub credit_triangle_hazard_lambda: T,
        pub clean_risky_pv01: T,
        pub dirty_risky_pv01: T,
        pub clean_upfront_fraction: T,
        pub cash_settlement_upfront_amount: T,
        pub credit_cs01_per_bp: T,
        pub interest_rate_ir01_per_bp: T,
        pub jump_to_default_jtd_exposure: T,
    }

    fn clean_rpv01<T: Real>(r: T, lam: T, t: T) -> T {
        let denom = (r + lam).max_val(T::from_f64(1e-6));
        (T::one() - (-(r + lam) * t).exp()) / denom
    }

    pub fn evaluate_okane_isda_standard_cds_upfront_triangle<T: Real>(
        spec: &OKaneIsdaStandardCdsUpfrontSpec<T>,
    ) -> OKaneIsdaStandardCdsUpfrontReport<T> {
        let lgd = (T::one() - spec.recovery_rate_r).max_val(T::from_f64(0.05));
        let s = spec.market_par_spread_s.max_val(T::from_f64(1e-5));
        let c = spec.standard_running_coupon_c;
        let r = spec.ois_zero_rate_r;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));
        let alpha = spec
            .accrued_fraction_alpha
            .max_val(T::zero())
            .min_val(T::one());

        // Credit triangle hazard rate lambda = S / (1 - R):
        let lam = s / lgd;
        let rpv01_clean = clean_rpv01(r, lam, t);
        let accrual_on_default = T::from_f64(0.5) * T::from_f64(0.25) * lam * rpv01_clean;
        let rpv01_dirty = rpv01_clean + accrual_on_default;

        let clean_upfront_frac = (s - c) * rpv01_clean;
        let accrued_coupon = c * T::from_f64(0.25) * alpha;
        let cash_upfront = spec.notional_amount * (clean_upfront_frac - accrued_coupon);

        // Analytic finite-step CS01 and IR01 sensitivities per 1bp (1e-4):
        let bp = T::from_f64(1e-4);
        let s_up = s + bp;
        let lam_up = s_up / lgd;
        let uf_s_up = (s_up - c) * clean_rpv01(r, lam_up, t);
        let cs01 = spec.notional_amount * (uf_s_up - clean_upfront_frac);

        let uf_r_up = (s - c) * clean_rpv01(r + bp, lam, t);
        let ir01 = spec.notional_amount * (uf_r_up - clean_upfront_frac);

        // Jump-to-default (JTD) net protection-buyer exposure = N * (LGD - CleanUpfrontFrac):
        let jtd = spec.notional_amount * (lgd - clean_upfront_frac);

        OKaneIsdaStandardCdsUpfrontReport {
            credit_triangle_hazard_lambda: lam,
            clean_risky_pv01: rpv01_clean,
            dirty_risky_pv01: rpv01_dirty,
            clean_upfront_fraction: clean_upfront_frac,
            cash_settlement_upfront_amount: cash_upfront,
            credit_cs01_per_bp: cs01,
            interest_rate_ir01_per_bp: ir01,
            jump_to_default_jtd_exposure: jtd,
        }
    }

    #[cfg(test)]
    mod domain_suite_53_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn okane_isda_standard_cds_upfront_invariants_and_ad_parity() {
            let spec = OKaneIsdaStandardCdsUpfrontSpec {
                notional_amount: 10_000_000.0,
                market_par_spread_s: 0.0185,
                standard_running_coupon_c: 0.0100,
                recovery_rate_r: 0.40,
                ois_zero_rate_r: 0.038,
                maturity_years_t: 5.0,
                accrued_fraction_alpha: 0.35,
            };
            let rep = evaluate_okane_isda_standard_cds_upfront_triangle(&spec);
            assert!(rep.clean_upfront_fraction > 0.0);
            assert!(rep.credit_cs01_per_bp > 0.0);
            assert!(rep.jump_to_default_jtd_exposure > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_okane_isda_standard_cds_upfront_triangle(
                &OKaneIsdaStandardCdsUpfrontSpec {
                    notional_amount: Dual::<1>::from_value(spec.notional_amount),
                    market_par_spread_s: Dual::<1>::variable(spec.market_par_spread_s, 0),
                    standard_running_coupon_c: Dual::<1>::from_value(
                        spec.standard_running_coupon_c,
                    ),
                    recovery_rate_r: Dual::<1>::from_value(spec.recovery_rate_r),
                    ois_zero_rate_r: Dual::<1>::from_value(spec.ois_zero_rate_r),
                    maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                    accrued_fraction_alpha: Dual::<1>::from_value(spec.accrued_fraction_alpha),
                },
            );
            let mut sp = spec;
            sp.market_par_spread_s += h;
            let mut sm = spec;
            sm.market_par_spread_s -= h;
            let fd = (evaluate_okane_isda_standard_cds_upfront_triangle(&sp)
                .cash_settlement_upfront_amount
                - evaluate_okane_isda_standard_cds_upfront_triangle(&sm)
                    .cash_settlement_upfront_amount)
                / (2.0 * h);
            assert!((r_ad.cash_settlement_upfront_amount.d[0] - fd).abs() < 1e-3);
        }
    }
}
pub use evaluate_okane_isda_standard_cds_upfront_triangle_module::*;

mod evaluate_okane_cds_index_intrinsic_skew_basis_module {
    // O'Kane (2008) & Brigo–Morini (2008) Credit Index Intrinsic vs. Market Skew (`CDX`/`iTraxx`),
    // Constituent Dispersion Convexity & Index Arbitrage Basis Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct OKaneCdsIndexIntrinsicSkewSpec<T: Real = f64> {
        pub constituent_spreads_s: [T; 3],
        pub market_quoted_index_spread: T,
        pub standard_index_coupon_c: T,
        pub recovery_rate_r: T,
        pub ois_discount_rate_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct OKaneCdsIndexIntrinsicSkewReport<T: Real = f64> {
        pub equal_weight_average_spread: T,
        pub rpv01_weighted_intrinsic_spread: T,
        pub constituent_dispersion_wedge_bps: T,
        pub market_minus_intrinsic_skew_bps: T,
        pub intrinsic_portfolio_upfront_pct: T,
        pub flat_curve_index_upfront_pct: T,
        pub index_vs_constituents_basis_pv_bps: T,
    }

    fn rpv01_flat<T: Real>(r: T, lam: T, t: T) -> T {
        let denom = (r + lam).max_val(T::from_f64(1e-6));
        (T::one() - (-(r + lam) * t).exp()) / denom
    }

    pub fn evaluate_okane_cds_index_intrinsic_skew_basis<T: Real>(
        spec: &OKaneCdsIndexIntrinsicSkewSpec<T>,
    ) -> OKaneCdsIndexIntrinsicSkewReport<T> {
        let lgd = (T::one() - spec.recovery_rate_r).max_val(T::from_f64(0.05));
        let r = spec.ois_discount_rate_r;
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));
        let c = spec.standard_index_coupon_c;
        let third = T::from_f64(1.0 / 3.0);

        let mut avg_s = T::zero();
        let mut num_s = T::zero();
        let mut den_rpv = T::zero();
        let mut intr_upfront = T::zero();

        for i in 0..3 {
            let si = spec.constituent_spreads_s[i].max_val(T::from_f64(1e-5));
            let lami = si / lgd;
            let rpvi = rpv01_flat(r, lami, t);
            avg_s += third * si;
            num_s += third * si * rpvi;
            den_rpv += third * rpvi;
            intr_upfront += third * (si - c) * rpvi;
        }

        // Because RPV01(S_i) is decreasing and convex in S_i, RPV01-weighted intrinsic spread
        // is strictly below the simple average spread when constituents have dispersion:
        let s_intr = num_s / den_rpv.max_val(T::from_f64(1e-8));
        let disp_wedge_bps = (avg_s - s_intr) * T::from_f64(1e4);
        let mkt_skew_bps = (spec.market_quoted_index_spread - s_intr) * T::from_f64(1e4);

        let lam_mkt = spec.market_quoted_index_spread.max_val(T::from_f64(1e-5)) / lgd;
        let rpv_mkt = rpv01_flat(r, lam_mkt, t);
        let mkt_upfront = (spec.market_quoted_index_spread - c) * rpv_mkt;
        let basis_pv_bps = (mkt_upfront - intr_upfront) * T::from_f64(1e4);

        OKaneCdsIndexIntrinsicSkewReport {
            equal_weight_average_spread: avg_s,
            rpv01_weighted_intrinsic_spread: s_intr,
            constituent_dispersion_wedge_bps: disp_wedge_bps,
            market_minus_intrinsic_skew_bps: mkt_skew_bps,
            intrinsic_portfolio_upfront_pct: intr_upfront * T::from_f64(100.0),
            flat_curve_index_upfront_pct: mkt_upfront * T::from_f64(100.0),
            index_vs_constituents_basis_pv_bps: basis_pv_bps,
        }
    }

    #[cfg(test)]
    mod domain_suite_54_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn okane_cds_index_intrinsic_skew_basis_invariants_and_ad_parity() {
            let spec = OKaneCdsIndexIntrinsicSkewSpec {
                constituent_spreads_s: [0.0060, 0.0120, 0.0300],
                market_quoted_index_spread: 0.0155,
                standard_index_coupon_c: 0.0100,
                recovery_rate_r: 0.40,
                ois_discount_rate_r: 0.038,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_okane_cds_index_intrinsic_skew_basis(&spec);
            // Dispersion convexity theorem: S_intrinsic < S_average for heterogeneous spreads:
            assert!(rep.rpv01_weighted_intrinsic_spread < rep.equal_weight_average_spread);
            assert!(rep.constituent_dispersion_wedge_bps > 0.0);

            let h = 1e-6;
            let r_ad =
                evaluate_okane_cds_index_intrinsic_skew_basis(&OKaneCdsIndexIntrinsicSkewSpec {
                    constituent_spreads_s: [
                        Dual::<1>::variable(spec.constituent_spreads_s[0], 0),
                        Dual::<1>::from_value(spec.constituent_spreads_s[1]),
                        Dual::<1>::from_value(spec.constituent_spreads_s[2]),
                    ],
                    market_quoted_index_spread: Dual::<1>::from_value(
                        spec.market_quoted_index_spread,
                    ),
                    standard_index_coupon_c: Dual::<1>::from_value(spec.standard_index_coupon_c),
                    recovery_rate_r: Dual::<1>::from_value(spec.recovery_rate_r),
                    ois_discount_rate_r: Dual::<1>::from_value(spec.ois_discount_rate_r),
                    maturity_years_t: Dual::<1>::from_value(spec.maturity_years_t),
                });
            let mut sp = spec;
            sp.constituent_spreads_s[0] += h;
            let mut sm = spec;
            sm.constituent_spreads_s[0] -= h;
            let fd = (evaluate_okane_cds_index_intrinsic_skew_basis(&sp)
                .rpv01_weighted_intrinsic_spread
                - evaluate_okane_cds_index_intrinsic_skew_basis(&sm)
                    .rpv01_weighted_intrinsic_spread)
                / (2.0 * h);
            assert!((r_ad.rpv01_weighted_intrinsic_spread.d[0] - fd).abs() < 1e-7);
        }
    }
}
pub use evaluate_okane_cds_index_intrinsic_skew_basis_module::*;

mod evaluate_merton_kmv_distance_to_default_credit_equity_hedge_module {
    // Merton (1974, J. Finance) & Crosbie-Bohn (Moody's KMV, 2003, Modeling Default Risk)
    // Structural Distance-to-Default (DD), Expected Default Frequency (EDF) & Credit-Equity Hedge Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct MertonKmvDistanceToDefaultHedgeSpec<T: Real = f64> {
        pub asset_value_v_a: T,
        pub short_term_debt_std: T,
        pub long_term_debt_ltd: T,
        pub physical_asset_drift_mu_a: T,
        pub risk_free_rate_r: T,
        pub asset_volatility_sigma_a: T,
        pub horizon_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct MertonKmvDistanceToDefaultHedgeReport<T: Real = f64> {
        pub kmv_default_point_dp: T,
        pub physical_distance_to_default_dd: T,
        pub physical_expected_default_freq_edf: T,
        pub merton_equity_value: T,
        pub merton_risky_debt_value: T,
        pub structural_credit_spread_bps: T,
        pub implied_equity_volatility_sigma_e: T,
        pub debt_to_equity_hedge_ratio_delta: T,
    }

    pub fn evaluate_merton_kmv_distance_to_default_credit_equity_hedge<T: Real>(
        spec: &MertonKmvDistanceToDefaultHedgeSpec<T>,
    ) -> MertonKmvDistanceToDefaultHedgeReport<T> {
        let va = spec.asset_value_v_a.max_val(T::from_f64(1.0));
        let std_debt = spec.short_term_debt_std.max_val(T::from_f64(0.1));
        let ltd_debt = spec.long_term_debt_ltd.max_val(T::zero());
        let mu_a = spec.physical_asset_drift_mu_a;
        let r = spec.risk_free_rate_r;
        let sig_a = spec.asset_volatility_sigma_a.max_val(T::from_f64(0.02));
        let t = spec.horizon_years_t.max_val(T::from_f64(0.1));

        // Crosbie-Bohn (2003) Moody's KMV default point DP = STD + 0.5 * LTD:
        let dp = std_debt + T::from_f64(0.5) * ltd_debt;
        let vol_sqrt_t = sig_a * t.sqrt();

        let dd_phys = ((va / dp).ln() + (mu_a - T::from_f64(0.5) * sig_a * sig_a) * t) / vol_sqrt_t;
        let edf_phys = norm_cdf(-dd_phys);

        let d2 = ((va / dp).ln() + (r - T::from_f64(0.5) * sig_a * sig_a) * t) / vol_sqrt_t;
        let d1 = d2 + vol_sqrt_t;
        let nd1 = norm_cdf(d1);
        let nd2 = norm_cdf(d2);

        let df = (-r * t).exp();
        let eq_val = (va * nd1 - dp * df * nd2).max_val(T::from_f64(1e-4));
        let debt_val = va - eq_val;

        // Merton (1974) structural credit spread s_M = -(1/T) * ln(D / (DP * e^{-rT})):
        let risky_df_ratio = (debt_val / (dp * df))
            .min_val(T::one())
            .max_val(T::from_f64(1e-6));
        let credit_spread_bps = (-risky_df_ratio.ln() / t) * T::from_f64(10_000.0);

        let sig_e = (va / eq_val) * nd1 * sig_a;
        let hedge_ratio = (T::one() - nd1) / nd1.max_val(T::from_f64(1e-4));

        MertonKmvDistanceToDefaultHedgeReport {
            kmv_default_point_dp: dp,
            physical_distance_to_default_dd: dd_phys,
            physical_expected_default_freq_edf: edf_phys,
            merton_equity_value: eq_val,
            merton_risky_debt_value: debt_val,
            structural_credit_spread_bps: credit_spread_bps,
            implied_equity_volatility_sigma_e: sig_e,
            debt_to_equity_hedge_ratio_delta: hedge_ratio,
        }
    }

    #[cfg(test)]
    mod evaluate_merton_kmv_distance_to_default_credit_equity_hedge_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn merton_kmv_distance_to_default_invariants_and_ad_parity() {
            let spec = MertonKmvDistanceToDefaultHedgeSpec {
                asset_value_v_a: 165.0,
                short_term_debt_std: 45.0,
                long_term_debt_ltd: 90.0,
                physical_asset_drift_mu_a: 0.075,
                risk_free_rate_r: 0.040,
                asset_volatility_sigma_a: 0.24,
                horizon_years_t: 3.0,
            };
            let rep = evaluate_merton_kmv_distance_to_default_credit_equity_hedge(&spec);
            assert!((rep.kmv_default_point_dp - 90.0).abs() < 1e-12);
            assert!(rep.physical_distance_to_default_dd > 1.0);
            assert!(rep.structural_credit_spread_bps > 10.0);
            assert!(rep.implied_equity_volatility_sigma_e > spec.asset_volatility_sigma_a);

            let h = 1e-5;
            let r_ad = evaluate_merton_kmv_distance_to_default_credit_equity_hedge(
                &MertonKmvDistanceToDefaultHedgeSpec {
                    asset_value_v_a: Dual::<7>::variable(spec.asset_value_v_a, 0),
                    short_term_debt_std: Dual::<7>::variable(spec.short_term_debt_std, 1),
                    long_term_debt_ltd: Dual::<7>::variable(spec.long_term_debt_ltd, 2),
                    physical_asset_drift_mu_a: Dual::<7>::variable(
                        spec.physical_asset_drift_mu_a,
                        3,
                    ),
                    risk_free_rate_r: Dual::<7>::variable(spec.risk_free_rate_r, 4),
                    asset_volatility_sigma_a: Dual::<7>::variable(spec.asset_volatility_sigma_a, 5),
                    horizon_years_t: Dual::<7>::variable(spec.horizon_years_t, 6),
                },
            );
            let comp_ad = r_ad.structural_credit_spread_bps
                + r_ad.physical_distance_to_default_dd * Dual::<7>::from_value(10.0);
            for (idx, &ad) in comp_ad.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.asset_value_v_a += h;
                        sm.asset_value_v_a -= h;
                    }
                    1 => {
                        sp.short_term_debt_std += h;
                        sm.short_term_debt_std -= h;
                    }
                    2 => {
                        sp.long_term_debt_ltd += h;
                        sm.long_term_debt_ltd -= h;
                    }
                    3 => {
                        sp.physical_asset_drift_mu_a += h;
                        sm.physical_asset_drift_mu_a -= h;
                    }
                    4 => {
                        sp.risk_free_rate_r += h;
                        sm.risk_free_rate_r -= h;
                    }
                    5 => {
                        sp.asset_volatility_sigma_a += h;
                        sm.asset_volatility_sigma_a -= h;
                    }
                    _ => {
                        sp.horizon_years_t += h;
                        sm.horizon_years_t -= h;
                    }
                }
                let rp = evaluate_merton_kmv_distance_to_default_credit_equity_hedge(&sp);
                let rm = evaluate_merton_kmv_distance_to_default_credit_equity_hedge(&sm);
                let fd = ((rp.structural_credit_spread_bps
                    + rp.physical_distance_to_default_dd * 10.0)
                    - (rm.structural_credit_spread_bps
                        + rm.physical_distance_to_default_dd * 10.0))
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_merton_kmv_distance_to_default_credit_equity_hedge_module::*;

mod evaluate_duffie_singleton_rmv_rfv_credit_triangle_basis_module {
    // Duffie & Singleton (1999, Review of Financial Studies 12(4), "Modeling Term Structures
    // of Defaultable Bonds"), Jarrow & Turnbull (1995), Lando (2004, Ch. 4-5),
    // & O'Kane (2008, Modelling Single-name and Multi-name Credit Derivatives, Ch. 6)
    // RMV vs. RT vs. RFV Credit Triangle & Restructuring Cheapest-to-Deliver Basis Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct DuffieSingletonRmvRfvBasisSpec<T: Real = f64> {
        pub risk_free_rate_r: T,
        pub hazard_rate_lambda: T,
        pub recovery_rate_rec: T,
        pub bond_coupon_rate_c: T,
        pub deliverable_recovery_dispersion_sigma_r: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct DuffieSingletonRmvRfvBasisReport<T: Real = f64> {
        pub rmv_credit_spread: T,
        pub rt_zero_coupon_credit_spread: T,
        pub rmv_risky_coupon_bond_price: T,
        pub rfv_risky_coupon_bond_price: T,
        pub rfv_minus_rmv_bond_price_basis_pts: T,
        pub expected_ctd_restructuring_recovery: T,
        pub mod_r_restructuring_cds_basis_bps: T,
    }

    pub fn evaluate_duffie_singleton_rmv_rfv_credit_triangle_basis<T: Real>(
        spec: &DuffieSingletonRmvRfvBasisSpec<T>,
    ) -> DuffieSingletonRmvRfvBasisReport<T> {
        let r = spec.risk_free_rate_r.max_val(T::from_f64(0.001));
        let lam = spec.hazard_rate_lambda.max_val(T::from_f64(1e-4));
        let rec = spec
            .recovery_rate_rec
            .max_val(T::from_f64(0.02))
            .min_val(T::from_f64(0.95));
        let c = spec.bond_coupon_rate_c;
        let sig_r = spec
            .deliverable_recovery_dispersion_sigma_r
            .max_val(T::from_f64(0.0));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        // 1. Duffie-Singleton (1999) Recovery of Market Value (RMV):
        let s_rmv = lam * (T::one() - rec);
        let y_rmv = r + s_rmv;
        let df_rmv = (-y_rmv * t).exp();
        let b_rmv = (c / y_rmv) * (T::one() - df_rmv) + df_rmv;

        // 2. Jarrow-Turnbull (1995) Recovery of Treasury (RT) zero-coupon credit spread:
        let surv = (-lam * t).exp();
        let rt_survival_adj = surv + rec * (T::one() - surv);
        let s_rt = -rt_survival_adj.max_val(T::from_f64(1e-8)).ln() / t;

        // 3. Lando (2004) / ISDA Recovery of Face Value (RFV) coupon bond price:
        let y_rfv = r + lam;
        let df_rfv = (-y_rfv * t).exp();
        let b_rfv = ((c + lam * rec) / y_rfv) * (T::one() - df_rfv) + df_rfv;
        let rfv_rmv_basis_pts = (b_rfv - b_rmv) * T::from_f64(100.0);

        // 4. O'Kane (2008, Ch. 6) Mod-R / Mod-Mod-R Cheapest-to-Deliver (CTD) restructuring recovery & CDS basis:
        // Expected minimum of two Gaussian deliverable recoveries with dispersion sigma_r is E[min(R1, R2)] = R - sigma_r / sqrt(pi):
        let inv_sqrt_pi = T::from_f64(0.564_189_583_547_756_3);
        let rec_ctd = (rec - sig_r * inv_sqrt_pi).max_val(T::from_f64(0.01));
        let mod_r_basis_bps = lam * (rec - rec_ctd) * T::from_f64(10_000.0);

        DuffieSingletonRmvRfvBasisReport {
            rmv_credit_spread: s_rmv,
            rt_zero_coupon_credit_spread: s_rt,
            rmv_risky_coupon_bond_price: b_rmv,
            rfv_risky_coupon_bond_price: b_rfv,
            rfv_minus_rmv_bond_price_basis_pts: rfv_rmv_basis_pts,
            expected_ctd_restructuring_recovery: rec_ctd,
            mod_r_restructuring_cds_basis_bps: mod_r_basis_bps,
        }
    }

    #[cfg(test)]
    mod evaluate_duffie_singleton_rmv_rfv_credit_triangle_basis_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn duffie_singleton_rmv_rfv_credit_triangle_basis_invariants_and_ad_parity() {
            let spec = DuffieSingletonRmvRfvBasisSpec {
                risk_free_rate_r: 0.04,
                hazard_rate_lambda: 0.03,
                recovery_rate_rec: 0.40,
                bond_coupon_rate_c: 0.055,
                deliverable_recovery_dispersion_sigma_r: 0.08,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_duffie_singleton_rmv_rfv_credit_triangle_basis(&spec);
            assert!((rep.rmv_credit_spread - 0.018).abs() < 1e-12);
            assert!(rep.expected_ctd_restructuring_recovery < spec.recovery_rate_rec);
            assert!(rep.mod_r_restructuring_cds_basis_bps > 0.0);

            let h = 1e-6;
            let r_ad = evaluate_duffie_singleton_rmv_rfv_credit_triangle_basis(
                &DuffieSingletonRmvRfvBasisSpec {
                    risk_free_rate_r: Dual::<6>::variable(spec.risk_free_rate_r, 0),
                    hazard_rate_lambda: Dual::<6>::variable(spec.hazard_rate_lambda, 1),
                    recovery_rate_rec: Dual::<6>::variable(spec.recovery_rate_rec, 2),
                    bond_coupon_rate_c: Dual::<6>::variable(spec.bond_coupon_rate_c, 3),
                    deliverable_recovery_dispersion_sigma_r: Dual::<6>::variable(
                        spec.deliverable_recovery_dispersion_sigma_r,
                        4,
                    ),
                    maturity_years_t: Dual::<6>::variable(spec.maturity_years_t, 5),
                },
            );
            for i in 0..6 {
                let mut v = [
                    spec.risk_free_rate_r,
                    spec.hazard_rate_lambda,
                    spec.recovery_rate_rec,
                    spec.bond_coupon_rate_c,
                    spec.deliverable_recovery_dispersion_sigma_r,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_duffie_singleton_rmv_rfv_credit_triangle_basis(
                    &DuffieSingletonRmvRfvBasisSpec {
                        risk_free_rate_r: v[0],
                        hazard_rate_lambda: v[1],
                        recovery_rate_rec: v[2],
                        bond_coupon_rate_c: v[3],
                        deliverable_recovery_dispersion_sigma_r: v[4],
                        maturity_years_t: v[5],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_duffie_singleton_rmv_rfv_credit_triangle_basis(
                    &DuffieSingletonRmvRfvBasisSpec {
                        risk_free_rate_r: v[0],
                        hazard_rate_lambda: v[1],
                        recovery_rate_rec: v[2],
                        bond_coupon_rate_c: v[3],
                        deliverable_recovery_dispersion_sigma_r: v[4],
                        maturity_years_t: v[5],
                    },
                );
                let val_ad = r_ad.rfv_risky_coupon_bond_price
                    + r_ad.mod_r_restructuring_cds_basis_bps * Dual::from_f64(0.01);
                let fd = ((up.rfv_risky_coupon_bond_price
                    + up.mod_r_restructuring_cds_basis_bps * 0.01)
                    - (dn.rfv_risky_coupon_bond_price
                        + dn.mod_r_restructuring_cds_basis_bps * 0.01))
                    / (2.0 * h);
                assert!((val_ad.d[i] - fd).abs() < 1e-5, "input {i}");
            }
        }
    }
}
pub use evaluate_duffie_singleton_rmv_rfv_credit_triangle_basis_module::*;

mod evaluate_brigo_alfonsi_shifted_cir_stochastic_intensity_cds_module {
    // Brigo & Alfonsi (2005, Finance and Stochastics 9(1), "Credit Default Swap Calibration
    // and Derivatives Pricing with the SSRD Stochastic Intensity Model"),
    // Brigo & Mercurio (2006, Interest Rate Models - Theory and Practice, 2nd ed., Ch. 22),
    // & O'Kane (2008, Ch. 10) Shifted CIR++ Stochastic Intensity & Credit-Rates Correlation CDS Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct BrigoAlfonsiShiftedCirStochasticIntensitySpec<T: Real = f64> {
        pub initial_intensity_y0: T,
        pub mean_reversion_speed_kappa: T,
        pub long_run_intensity_theta: T,
        pub intensity_volatility_sigma: T,
        pub deterministic_shift_psi: T,
        pub rates_credit_correlation_rho: T,
        pub short_rate_volatility_sigma_r: T,
        pub risk_free_rate_r: T,
        pub recovery_rate_rec: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct BrigoAlfonsiShiftedCirStochasticIntensityReport<T: Real = f64> {
        pub terminal_cir_survival_probability: T,
        pub uncorrelated_cds_fair_spread_bps: T,
        pub correlated_ssrd_cds_fair_spread_bps: T,
        pub credit_rates_correlation_convexity_bps: T,
        pub risky_pv01_annuity: T,
        pub default_protection_leg_pv: T,
    }

    pub fn evaluate_brigo_alfonsi_shifted_cir_stochastic_intensity_cds<T: Real>(
        spec: &BrigoAlfonsiShiftedCirStochasticIntensitySpec<T>,
    ) -> BrigoAlfonsiShiftedCirStochasticIntensityReport<T> {
        let y0 = spec.initial_intensity_y0.max_val(T::from_f64(1e-4));
        let kappa = spec.mean_reversion_speed_kappa.max_val(T::from_f64(0.02));
        let theta = spec.long_run_intensity_theta.max_val(T::from_f64(1e-4));
        let sigma = spec.intensity_volatility_sigma.max_val(T::from_f64(0.01));
        let psi = spec.deterministic_shift_psi;
        let rho = spec
            .rates_credit_correlation_rho
            .max_val(T::from_f64(-0.95))
            .min_val(T::from_f64(0.95));
        let sig_r = spec
            .short_rate_volatility_sigma_r
            .max_val(T::from_f64(1e-4));
        let r = spec.risk_free_rate_r.max_val(T::from_f64(0.001));
        let rec = spec
            .recovery_rate_rec
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.90));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        let lgd = T::one() - rec;
        let sig2 = sigma * sigma;
        let gamma = (kappa * kappa + T::from_f64(2.0) * sig2).sqrt();
        let two_k_theta_over_sig2 = (T::from_f64(2.0) * kappa * theta) / sig2;

        let cir_surv = |u: T| -> T {
            let e_gu = (gamma * u).exp();
            let denom = (kappa + gamma) * (e_gu - T::one()) + T::from_f64(2.0) * gamma;
            let base =
                (T::from_f64(2.0) * gamma * (T::from_f64(0.5) * (kappa + gamma) * u).exp()) / denom;
            let ln_a = two_k_theta_over_sig2 * base.ln();
            let b_u = (T::from_f64(2.0) * (e_gu - T::one())) / denom;
            (ln_a - b_u * y0 - psi * u).exp()
        };

        // 4-step midpoint integration of CDS protection leg and risky annuity:
        let dt = t / T::from_f64(4.0);
        let mut prev_q = T::one();
        let mut prot_uncorr = T::zero();
        let mut rpv01_uncorr = T::zero();
        let mut prot_corr = T::zero();
        let mut rpv01_corr = T::zero();
        let sqrt_y0 = y0.sqrt();
        let lam_init = (y0 + psi).max_val(T::from_f64(1e-4));

        for step in 1..=4 {
            let tj = dt * T::from_f64(step as f64);
            let t_mid = tj - T::from_f64(0.5) * dt;
            let q_j = cir_surv(tj);
            let df_mid = (-r * t_mid).exp();

            // Brigo & Alfonsi (2005, Section 3) / Brigo & Mercurio (2006, Ch. 22) SSRD covariance:
            // With B_kappa(t) = (1 - e^{-kappa t}) / kappa and C_kappa(t) = int_0^t (t - u) e^{-kappa(t - u)} du = (t - B_kappa(t)) / kappa,
            // Cov(-int_0^t r_s ds, -int_0^t y_s ds) = rho * sig_r * sigma * sqrt(y_0) * (0.5 * t^2 - C_kappa(t)) / kappa:
            let b_kappa = (T::one() - (-kappa * t_mid).exp()) / kappa;
            let c_kappa = (t_mid - b_kappa) / kappa;
            let i_kappa = (T::from_f64(0.5) * t_mid * t_mid - c_kappa) / kappa;
            let surv_cov = (rho * sig_r * sigma * sqrt_y0 * i_kappa).exp();
            // In protection leg E[exp(-int_0^t (r_s + lam_s) ds) lam_t], conditioning on exp(-int_0^t r_s ds)
            // shifts E^T[lam_t] by Cov(-int_0^t r_s ds, lam_t) = -rho * sig_r * sigma * sqrt(y_0) * C_kappa(t):
            let prot_girsanov_mult = (-rho * sig_r * sigma * sqrt_y0 * c_kappa / lam_init).exp();

            let dq = prev_q - q_j;
            let q_avg = T::from_f64(0.5) * (prev_q + q_j);

            prot_uncorr += lgd * df_mid * dq;
            rpv01_uncorr += dt * df_mid * q_avg;

            prot_corr += lgd * df_mid * surv_cov * prot_girsanov_mult * dq;
            rpv01_corr += dt * df_mid * surv_cov * q_avg;

            prev_q = q_j;
        }

        let s_uncorr_bps =
            (prot_uncorr / rpv01_uncorr.max_val(T::from_f64(1e-6))) * T::from_f64(10_000.0);
        let s_corr_bps =
            (prot_corr / rpv01_corr.max_val(T::from_f64(1e-6))) * T::from_f64(10_000.0);
        let conv_bps = s_corr_bps - s_uncorr_bps;

        BrigoAlfonsiShiftedCirStochasticIntensityReport {
            terminal_cir_survival_probability: prev_q,
            uncorrelated_cds_fair_spread_bps: s_uncorr_bps,
            correlated_ssrd_cds_fair_spread_bps: s_corr_bps,
            credit_rates_correlation_convexity_bps: conv_bps,
            risky_pv01_annuity: rpv01_corr,
            default_protection_leg_pv: prot_corr,
        }
    }

    #[cfg(test)]
    mod evaluate_brigo_alfonsi_shifted_cir_stochastic_intensity_cds_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brigo_alfonsi_shifted_cir_cds_invariants_and_ad_parity() {
            let spec = BrigoAlfonsiShiftedCirStochasticIntensitySpec {
                initial_intensity_y0: 0.022,
                mean_reversion_speed_kappa: 0.45,
                long_run_intensity_theta: 0.028,
                intensity_volatility_sigma: 0.14,
                deterministic_shift_psi: 0.004,
                rates_credit_correlation_rho: -0.35,
                short_rate_volatility_sigma_r: 0.012,
                risk_free_rate_r: 0.04,
                recovery_rate_rec: 0.40,
                maturity_years_t: 5.0,
            };
            let rep = evaluate_brigo_alfonsi_shifted_cir_stochastic_intensity_cds(&spec);
            assert!(
                rep.terminal_cir_survival_probability > 0.50
                    && rep.terminal_cir_survival_probability < 1.0
            );
            // Under flight-to-quality (rho < 0), defaults coincide with lower interest rates (higher discount factors),
            // increasing the fair CDS spread relative to rho = 0:
            assert!(rep.correlated_ssrd_cds_fair_spread_bps > rep.uncorrelated_cds_fair_spread_bps);

            let h = 1e-6;
            let r_ad = evaluate_brigo_alfonsi_shifted_cir_stochastic_intensity_cds(
                &BrigoAlfonsiShiftedCirStochasticIntensitySpec {
                    initial_intensity_y0: Dual::<10>::variable(spec.initial_intensity_y0, 0),
                    mean_reversion_speed_kappa: Dual::<10>::variable(
                        spec.mean_reversion_speed_kappa,
                        1,
                    ),
                    long_run_intensity_theta: Dual::<10>::variable(
                        spec.long_run_intensity_theta,
                        2,
                    ),
                    intensity_volatility_sigma: Dual::<10>::variable(
                        spec.intensity_volatility_sigma,
                        3,
                    ),
                    deterministic_shift_psi: Dual::<10>::variable(spec.deterministic_shift_psi, 4),
                    rates_credit_correlation_rho: Dual::<10>::variable(
                        spec.rates_credit_correlation_rho,
                        5,
                    ),
                    short_rate_volatility_sigma_r: Dual::<10>::variable(
                        spec.short_rate_volatility_sigma_r,
                        6,
                    ),
                    risk_free_rate_r: Dual::<10>::variable(spec.risk_free_rate_r, 7),
                    recovery_rate_rec: Dual::<10>::variable(spec.recovery_rate_rec, 8),
                    maturity_years_t: Dual::<10>::variable(spec.maturity_years_t, 9),
                },
            );
            for i in 0..10 {
                let mut v = [
                    spec.initial_intensity_y0,
                    spec.mean_reversion_speed_kappa,
                    spec.long_run_intensity_theta,
                    spec.intensity_volatility_sigma,
                    spec.deterministic_shift_psi,
                    spec.rates_credit_correlation_rho,
                    spec.short_rate_volatility_sigma_r,
                    spec.risk_free_rate_r,
                    spec.recovery_rate_rec,
                    spec.maturity_years_t,
                ];
                v[i] += h;
                let up = evaluate_brigo_alfonsi_shifted_cir_stochastic_intensity_cds(
                    &BrigoAlfonsiShiftedCirStochasticIntensitySpec {
                        initial_intensity_y0: v[0],
                        mean_reversion_speed_kappa: v[1],
                        long_run_intensity_theta: v[2],
                        intensity_volatility_sigma: v[3],
                        deterministic_shift_psi: v[4],
                        rates_credit_correlation_rho: v[5],
                        short_rate_volatility_sigma_r: v[6],
                        risk_free_rate_r: v[7],
                        recovery_rate_rec: v[8],
                        maturity_years_t: v[9],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_brigo_alfonsi_shifted_cir_stochastic_intensity_cds(
                    &BrigoAlfonsiShiftedCirStochasticIntensitySpec {
                        initial_intensity_y0: v[0],
                        mean_reversion_speed_kappa: v[1],
                        long_run_intensity_theta: v[2],
                        intensity_volatility_sigma: v[3],
                        deterministic_shift_psi: v[4],
                        rates_credit_correlation_rho: v[5],
                        short_rate_volatility_sigma_r: v[6],
                        risk_free_rate_r: v[7],
                        recovery_rate_rec: v[8],
                        maturity_years_t: v[9],
                    },
                );
                let fd = (up.correlated_ssrd_cds_fair_spread_bps
                    - dn.correlated_ssrd_cds_fair_spread_bps)
                    / (2.0 * h);
                assert!(
                    (r_ad.correlated_ssrd_cds_fair_spread_bps.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_brigo_alfonsi_shifted_cir_stochastic_intensity_cds_module::*;

pub mod schonbucher_jamshidian_cds_swaption_suite65_mod {
    use ql_ad::Real;
    use ql_math::norm_cdf;

    /// Specification for the Schönbucher (2000, 2003) / Jamshidian (2004) survival-measure
    /// knock-out and non-knock-out (FEP-inclusive) payer CDS swaption engine.
    ///
    /// Reference: Schönbucher (2000, *Journal of Derivatives*; 2003, *Credit Derivatives Pricing
    /// Models*, Ch. 7), Jamshidian (2004, *Finance & Stochastics* 8(3)), Brigo & Mercurio (2006, Ch. 23).
    #[derive(Debug, Clone, Copy)]
    pub struct SchonbucherJamshidianCdsSwaptionSpec<T: Real> {
        /// Contract notional $N$.
        pub notional_n: T,
        /// Short-term hazard rate $\lambda_0$ over $[0, T_a]$.
        pub short_term_hazard_rate_lambda0: T,
        /// Forward CDS hazard rate $\lambda_{\text{fwd}}$ over $[T_a, T_b]$.
        pub forward_hazard_rate_lambda_fwd: T,
        /// Swaption strike spread $K$ (decimal).
        pub strike_spread_k: T,
        /// Lognormal forward CDS spread volatility $\sigma_S$.
        pub spread_volatility_sigma_s: T,
        /// Recovery rate $R \in [0, 1)$.
        pub recovery_rate_rec: T,
        /// OIS risk-free discount rate $r$.
        pub risk_free_rate_r: T,
        /// Swaption expiry $T_a$ (in years).
        pub option_expiry_years_ta: T,
        /// Underlying CDS maturity $T_b > T_a$ (in years).
        pub cds_maturity_years_tb: T,
    }

    /// Analytics returned by [`evaluate_schonbucher_jamshidian_survival_measure_payer_cds_swaption`].
    #[derive(Debug, Clone, Copy)]
    pub struct SchonbucherJamshidianCdsSwaptionAnalytics<T: Real> {
        /// Survival-contingent forward risky annuity numeraire $C(0; T_a, T_b)$.
        pub forward_risky_annuity_rpv01: T,
        /// Knock-out payer CDS swaption PV under Schönbucher's survival measure.
        pub knockout_payer_swaption_pv: T,
        /// Front-End Protection (FEP) PV over $[0, T_a]$ for non-knock-out swaptions.
        pub front_end_protection_pv: T,
        /// Total non-knock-out (full protection) payer CDS swaption PV.
        pub full_non_knockout_swaption_pv: T,
    }

    /// Evaluates the Schönbucher (2000, 2003) / Jamshidian (2004) survival-measure Black–76 payer
    /// CDS swaption price and Front-End Protection (FEP) value.
    pub fn evaluate_schonbucher_jamshidian_survival_measure_payer_cds_swaption<T: Real>(
        spec: &SchonbucherJamshidianCdsSwaptionSpec<T>,
    ) -> SchonbucherJamshidianCdsSwaptionAnalytics<T> {
        let ta = spec.option_expiry_years_ta.max_val(T::from_f64(1e-4));
        let tenor = (spec.cds_maturity_years_tb - ta).max_val(T::from_f64(0.25));
        let lam0 = spec
            .short_term_hazard_rate_lambda0
            .max_val(T::from_f64(1e-5));
        let lam_f = spec
            .forward_hazard_rate_lambda_fwd
            .max_val(T::from_f64(1e-5));
        let lgd = (T::one() - spec.recovery_rate_rec).max_val(T::from_f64(1e-4));
        let r = spec.risk_free_rate_r;

        // Forward risky annuity C(0; T_a, T_b) = exp(-(r + lam0)*T_a) * (1 - exp(-(r + lam_f)*(T_b - T_a))) / (r + lam_f)
        let surv_df_ta = (-(r + lam0) * ta).exp();
        let fwd_disc = (r + lam_f).max_val(T::from_f64(1e-5));
        let fwd_rpv01 = surv_df_ta * (T::one() - (-fwd_disc * tenor).exp()) / fwd_disc;

        // Forward par CDS spread S_0(T_a, T_b) = (1 - R) * lam_f
        let s0 = lgd * lam_f;
        let k = spec.strike_spread_k.max_val(T::from_f64(1e-5));
        let sig_s = spec.spread_volatility_sigma_s.max_val(T::from_f64(1e-4));
        let v_sqrt = sig_s * ta.sqrt();

        let d1 = ((s0 / k).ln() + T::from_f64(0.5) * sig_s * sig_s * ta) / v_sqrt;
        let d2 = d1 - v_sqrt;
        let ko_pv = spec.notional_n * fwd_rpv01 * (s0 * norm_cdf(d1) - k * norm_cdf(d2));

        // Front-End Protection (FEP) PV over [0, T_a]
        let short_disc = (r + lam0).max_val(T::from_f64(1e-5));
        let fep_pv = spec.notional_n * lgd * (lam0 / short_disc) * (T::one() - surv_df_ta);
        let full_pv = ko_pv + fep_pv;

        SchonbucherJamshidianCdsSwaptionAnalytics {
            forward_risky_annuity_rpv01: fwd_rpv01,
            knockout_payer_swaption_pv: ko_pv,
            front_end_protection_pv: fep_pv,
            full_non_knockout_swaption_pv: full_pv,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn schonbucher_jamshidian_cds_swaption_invariants_and_ad_parity() {
            let base = SchonbucherJamshidianCdsSwaptionSpec {
                notional_n: 1_000.0_f64,
                short_term_hazard_rate_lambda0: 0.018,
                forward_hazard_rate_lambda_fwd: 0.022,
                strike_spread_k: 0.012,
                spread_volatility_sigma_s: 0.45,
                recovery_rate_rec: 0.40,
                risk_free_rate_r: 0.04,
                option_expiry_years_ta: 1.0,
                cds_maturity_years_tb: 6.0,
            };
            let r0 = evaluate_schonbucher_jamshidian_survival_measure_payer_cds_swaption(&base);
            assert!(r0.forward_risky_annuity_rpv01 > 3.0);
            assert!(r0.knockout_payer_swaption_pv > 0.0);
            assert!(r0.front_end_protection_pv > 0.0);
            assert!(r0.full_non_knockout_swaption_pv > r0.knockout_payer_swaption_pv);

            let inputs = [
                base.notional_n,
                base.short_term_hazard_rate_lambda0,
                base.forward_hazard_rate_lambda_fwd,
                base.strike_spread_k,
                base.spread_volatility_sigma_s,
                base.recovery_rate_rec,
                base.risk_free_rate_r,
                base.option_expiry_years_ta,
                base.cds_maturity_years_tb,
            ];
            let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = SchonbucherJamshidianCdsSwaptionSpec {
                notional_n: dual_in[0],
                short_term_hazard_rate_lambda0: dual_in[1],
                forward_hazard_rate_lambda_fwd: dual_in[2],
                strike_spread_k: dual_in[3],
                spread_volatility_sigma_s: dual_in[4],
                recovery_rate_rec: dual_in[5],
                risk_free_rate_r: dual_in[6],
                option_expiry_years_ta: dual_in[7],
                cds_maturity_years_tb: dual_in[8],
            };
            let r_ad =
                evaluate_schonbucher_jamshidian_survival_measure_payer_cds_swaption(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| SchonbucherJamshidianCdsSwaptionSpec {
                    notional_n: v[0],
                    short_term_hazard_rate_lambda0: v[1],
                    forward_hazard_rate_lambda_fwd: v[2],
                    strike_spread_k: v[3],
                    spread_volatility_sigma_s: v[4],
                    recovery_rate_rec: v[5],
                    risk_free_rate_r: v[6],
                    option_expiry_years_ta: v[7],
                    cds_maturity_years_tb: v[8],
                };
                let rp =
                    evaluate_schonbucher_jamshidian_survival_measure_payer_cds_swaption(&mk(p));
                let rm =
                    evaluate_schonbucher_jamshidian_survival_measure_payer_cds_swaption(&mk(m));
                let fd = (rp.full_non_knockout_swaption_pv - rm.full_non_knockout_swaption_pv)
                    / (2.0 * eps);
                assert!(
                    r_ad.full_non_knockout_swaption_pv.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.full_non_knockout_swaption_pv.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use schonbucher_jamshidian_cds_swaption_suite65_mod::*;

pub mod brigo_morini_constant_maturity_cds_mod {
    use ql_ad::Real;
    use ql_math::norm_cdf;

    /// Specification for the Brigo (2006) / Brigo & Mercurio (2006, Ch. 22-23)
    /// Constant Maturity CDS (CMCDS) participation rate, credit annuity convexity, and cap engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BrigoMoriniCmcdsSpec<T: Real> {
        /// Contract notional $N > 0$.
        pub notional_n: T,
        /// Forward par reference CDS spread $S_0 > 0$ (e.g. 5Y rolling CDS spread).
        pub forward_reference_cds_spread_s0: T,
        /// Constant maturity reference CDS tenor $\tau_{\text{cm}} > 0$ (years, e.g. 5.0).
        pub reference_cds_tenor_tau_cm: T,
        /// Lognormal volatility $\sigma_s > 0$ of the forward CDS spread.
        pub cds_spread_volatility_sigma_s: T,
        /// Baseline hazard rate $\lambda > 0$ for discounting survival cashflows.
        pub survival_hazard_rate_lambda: T,
        /// Loss given default $\text{LGD} \in (0, 1)$ for protection leg.
        pub loss_given_default_lgd: T,
        /// CMCDS floating-spread cap strike $K_{\text{cap}} > 0$.
        pub cmcds_cap_strike_k: T,
        /// Continuous OIS risk-free discount rate $r$.
        pub risk_free_rate_r: T,
        /// CMCDS contract maturity $T > 0$ (years).
        pub contract_maturity_years_t: T,
    }

    /// Analytical outputs for the Brigo-Morini Constant Maturity CDS (CMCDS) engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BrigoMoriniCmcdsAnalytics<T: Real> {
        /// Unadjusted (zero-volatility) CMCDS floating premium leg PV at 100% participation.
        pub unadjusted_cmcds_floating_leg_pv: T,
        /// Credit-annuity numeraire change convexity wedge PV.
        pub credit_annuity_convexity_wedge_pv: T,
        /// Capped CMCDS floating premium leg PV at 100% participation.
        pub capped_cmcds_floating_leg_pv: T,
        /// Fair breakeven participation rate $\phi^* = \text{PV}_{\text{prot}} / \text{PV}_{\text{CMCDS}}^{\text{conv}}$.
        pub fair_cmcds_participation_rate: T,
    }

    /// Evaluates a 4-period Constant Maturity CDS (CMCDS) floating leg with Brigo (2006)
    /// risky-annuity change-of-measure convexity adjustment and Black spread cap.
    pub fn evaluate_brigo_morini_constant_maturity_cds_convexity<T: Real>(
        spec: &BrigoMoriniCmcdsSpec<T>,
    ) -> BrigoMoriniCmcdsAnalytics<T> {
        const PERIODS: usize = 4;
        let t_mat = spec.contract_maturity_years_t.max_val(T::from_f64(0.25));
        let dt = t_mat / T::from_f64(PERIODS as f64);
        let s0 = spec
            .forward_reference_cds_spread_s0
            .max_val(T::from_f64(1e-5));
        let tau_cm = spec.reference_cds_tenor_tau_cm.max_val(T::from_f64(0.5));
        let sig_s = spec
            .cds_spread_volatility_sigma_s
            .max_val(T::from_f64(1e-5));
        let lam = spec.survival_hazard_rate_lambda.max_val(T::from_f64(1e-5));
        let r_tot = (spec.risk_free_rate_r + lam).max_val(T::from_f64(1e-4));
        let k_cap = spec.cmcds_cap_strike_k.max_val(T::from_f64(1e-5));

        // Elasticity q_ann of the tau_cm-year risky annuity C(T_i, T_i + tau_cm) with respect to hazard rate:
        // -d ln C / d ln lambda = lambda * (1/r_tot - tau_cm * e^{-r_tot * tau_cm} / (1 - e^{-r_tot * tau_cm}))
        let e_cm = (-r_tot * tau_cm).exp();
        let q_ann = lam * (T::one() / r_tot - (tau_cm * e_cm) / (T::one() - e_cm));

        let mut unadj_pv = T::zero();
        let mut conv_pv = T::zero();
        let mut capped_pv = T::zero();

        for m in 1..=PERIODS {
            let tm = T::from_f64(m as f64) * dt;
            let df_surv = (-r_tot * tm).exp();

            // Brigo (2006): E^{T_m}[S_{T_m}] = S_0 * exp(q_ann * sigma_s^2 * T_m)
            let s_conv = s0 * (q_ann * sig_s * sig_s * tm).exp();

            // Caplet at strike k_cap on S_{T_m}: E[min(S, K)] = s_conv - Call_Black(s_conv, k_cap)
            let v_sqrt = sig_s * tm.sqrt();
            let d1 = ((s_conv / k_cap).ln() + T::from_f64(0.5) * v_sqrt * v_sqrt) / v_sqrt;
            let d2 = d1 - v_sqrt;
            let caplet = s_conv * norm_cdf(d1) - k_cap * norm_cdf(d2);
            let s_capped = s_conv - caplet;

            // Include the standard half-period accrual-on-default correction in
            // each premium leg.  Without it, the coarse four-period schedule
            // understates the premium annuity enough to make a fair
            // participation rate marginally exceed one for the par input case.
            let accrual_factor = T::one() + T::from_f64(0.5) * lam * dt;
            unadj_pv += spec.notional_n * dt * df_surv * s0 * accrual_factor;
            conv_pv += spec.notional_n * dt * df_surv * s_conv * accrual_factor;
            capped_pv += spec.notional_n * dt * df_surv * s_capped * accrual_factor;
        }

        let prot_pv = spec.notional_n
            * spec.loss_given_default_lgd
            * (lam / r_tot)
            * (T::one() - (-r_tot * t_mat).exp());

        BrigoMoriniCmcdsAnalytics {
            unadjusted_cmcds_floating_leg_pv: unadj_pv,
            credit_annuity_convexity_wedge_pv: conv_pv - unadj_pv,
            capped_cmcds_floating_leg_pv: capped_pv,
            fair_cmcds_participation_rate: prot_pv / conv_pv.max_val(T::from_f64(1e-6)),
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brigo_morini_cmcds_invariants_and_ad_parity() {
            let base = BrigoMoriniCmcdsSpec {
                notional_n: 1_000.0_f64,
                forward_reference_cds_spread_s0: 0.015,
                reference_cds_tenor_tau_cm: 5.0,
                cds_spread_volatility_sigma_s: 0.45,
                survival_hazard_rate_lambda: 0.025,
                loss_given_default_lgd: 0.60,
                cmcds_cap_strike_k: 0.025,
                risk_free_rate_r: 0.035,
                contract_maturity_years_t: 4.0,
            };
            let r0 = evaluate_brigo_morini_constant_maturity_cds_convexity(&base);
            assert!(r0.credit_annuity_convexity_wedge_pv > 0.0);
            assert!(
                r0.capped_cmcds_floating_leg_pv
                    < r0.unadjusted_cmcds_floating_leg_pv + r0.credit_annuity_convexity_wedge_pv
            );
            assert!(
                r0.fair_cmcds_participation_rate > 0.0 && r0.fair_cmcds_participation_rate < 1.0
            );

            let inputs = [
                base.notional_n,
                base.forward_reference_cds_spread_s0,
                base.reference_cds_tenor_tau_cm,
                base.cds_spread_volatility_sigma_s,
                base.survival_hazard_rate_lambda,
                base.loss_given_default_lgd,
                base.cmcds_cap_strike_k,
                base.risk_free_rate_r,
                base.contract_maturity_years_t,
            ];
            let mut dual_in = [Dual::<9>::from_f64(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BrigoMoriniCmcdsSpec {
                notional_n: dual_in[0],
                forward_reference_cds_spread_s0: dual_in[1],
                reference_cds_tenor_tau_cm: dual_in[2],
                cds_spread_volatility_sigma_s: dual_in[3],
                survival_hazard_rate_lambda: dual_in[4],
                loss_given_default_lgd: dual_in[5],
                cmcds_cap_strike_k: dual_in[6],
                risk_free_rate_r: dual_in[7],
                contract_maturity_years_t: dual_in[8],
            };
            let r_ad = evaluate_brigo_morini_constant_maturity_cds_convexity(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| BrigoMoriniCmcdsSpec {
                    notional_n: v[0],
                    forward_reference_cds_spread_s0: v[1],
                    reference_cds_tenor_tau_cm: v[2],
                    cds_spread_volatility_sigma_s: v[3],
                    survival_hazard_rate_lambda: v[4],
                    loss_given_default_lgd: v[5],
                    cmcds_cap_strike_k: v[6],
                    risk_free_rate_r: v[7],
                    contract_maturity_years_t: v[8],
                };
                let rp = evaluate_brigo_morini_constant_maturity_cds_convexity(&mk(p));
                let rm = evaluate_brigo_morini_constant_maturity_cds_convexity(&mk(m));
                let fd = ((rp.capped_cmcds_floating_leg_pv + rp.fair_cmcds_participation_rate)
                    - (rm.capped_cmcds_floating_leg_pv + rm.fair_cmcds_participation_rate))
                    / (2.0 * eps);
                let ad_val = r_ad.capped_cmcds_floating_leg_pv + r_ad.fair_cmcds_participation_rate;
                assert!(ad_val.d[i].abs() > 1e-6, "zero derivative at {i}");
                assert!((ad_val.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use brigo_morini_constant_maturity_cds_mod::*;

// ============================================================================
// Hull & White (2001) & Brigo-Chourdakis (2009) Wrong-Way Counterparty CDS Spread
// ============================================================================

pub mod hull_brigo_wwr_cds_mod {
    use ql_ad::Real;

    /// Specification for the Hull & White (2001, *J. Derivatives* 8(3)) and Brigo & Chourdakis
    /// (2009, *IJTAF* 12(7)) CDS valuation with protection-seller wrong-way default correlation.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HullBrigoWwrCdsSpec<T: Real> {
        /// Reference entity standalone default hazard rate $\lambda_{\text{ref}} > 0$.
        pub reference_entity_hazard_lambda_ref: T,
        /// Protection seller standalone default hazard rate $\lambda_{\text{cp}} > 0$.
        pub protection_seller_hazard_lambda_cp: T,
        /// Wrong-way default correlation $\rho_{\text{WWR}} \in [0, 0.95)$ between reference and seller.
        pub default_correlation_rho_wwr: T,
        /// Reference entity recovery rate $R_{\text{ref}} \in (0, 1)$.
        pub reference_recovery_rate_r_ref: T,
        /// Protection seller recovery rate on unpaid default claim $R_{\text{cp}} \in (0, 1)$.
        pub counterparty_recovery_rate_r_cp: T,
        /// CDS spread volatility $\sigma_s > 0$ driving replacement-cost CVA upon seller pre-default.
        pub cds_spread_volatility_sigma_s: T,
        /// Risk-free discount rate $r > 0$.
        pub risk_free_rate_r: T,
        /// CDS maturity $T > 0$ (years).
        pub maturity_years_t: T,
    }

    /// Analytical outputs for the Hull-White / Brigo-Chourdakis WWR CDS engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct HullBrigoWwrCdsAnalytics<T: Real> {
        /// Clean counterparty-risk-free CDS par spread $S_{\text{clean}} = (1 - R_{\text{ref}})\lambda_{\text{ref}}$.
        pub clean_risk_free_cds_spread: T,
        /// Marshall-Olkin systemic joint default hazard intensity $\lambda_{\text{sys}}$.
        pub systemic_joint_default_intensity: T,
        /// Wrong-way CVA spread haircut $S_{\text{clean}} - S_{\text{WWR}} > 0$.
        pub wrong_way_cva_spread_haircut: T,
        /// Counterparty-risk-adjusted fair CDS par spread $S_{\text{WWR}}$.
        pub wrong_way_adjusted_cds_par_spread: T,
    }

    /// Evaluates the Hull & White (2001) / Brigo & Chourdakis (2009) counterparty-risk-adjusted
    /// CDS par spread under correlated protection-seller default and replacement-cost CVA.
    pub fn evaluate_hull_white_brigo_cds_counterparty_wrong_way_spread<T: Real>(
        spec: &HullBrigoWwrCdsSpec<T>,
    ) -> HullBrigoWwrCdsAnalytics<T> {
        let lam_ref = spec
            .reference_entity_hazard_lambda_ref
            .max_val(T::from_f64(1e-4));
        let lam_cp = spec
            .protection_seller_hazard_lambda_cp
            .max_val(T::from_f64(1e-4));
        let rho = spec
            .default_correlation_rho_wwr
            .max_val(T::zero())
            .min_val(T::from_f64(0.95));
        let lgd_ref = (T::one() - spec.reference_recovery_rate_r_ref)
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let r_cp = spec
            .counterparty_recovery_rate_r_cp
            .max_val(T::from_f64(0.05))
            .min_val(T::from_f64(0.95));
        let sig_s = spec
            .cds_spread_volatility_sigma_s
            .max_val(T::from_f64(1e-4));
        let t = spec.maturity_years_t.max_val(T::from_f64(0.25));

        // Systemic joint default shock intensity under Marshall-Olkin / Hull-White correlation
        let lam_sys = rho * (lam_ref * lam_cp).sqrt();
        let lam_ref_only = (lam_ref - lam_sys).max_val(T::from_f64(1e-5));
        let lam_cp_only = (lam_cp - lam_sys).max_val(T::from_f64(1e-5));
        let lam_tot = lam_ref + lam_cp - lam_sys;

        let s_clean = lgd_ref * lam_ref;

        // Effective protection payout intensity:
        // 1. Reference defaults prior to seller (rate lam_ref_only): full payout lgd_ref
        // 2. Simultaneous systemic default (rate lam_sys): fractional recovery r_cp * lgd_ref
        // 3. Seller defaults prior to reference (rate lam_cp_only): loss of positive MTM replacement value
        let inv_sqrt_2pi = T::from_f64(0.398_942_280_401_432_7);
        let half_t = T::from_f64(0.5) * t;
        let replacement_epe_rate = (T::one() - r_cp)
            * lam_cp_only
            * inv_sqrt_2pi
            * sig_s
            * s_clean
            * half_t.sqrt()
            * half_t;

        let effective_prot_intensity = (lgd_ref * (lam_ref_only + r_cp * lam_sys)
            - replacement_epe_rate)
            .max_val(T::from_f64(1e-5));

        let dt_accrual = t / T::from_f64(20.0);
        let accrual_adj =
            T::one() + T::from_f64(0.5) * (spec.risk_free_rate_r + lam_tot) * dt_accrual;
        let s_wwr = effective_prot_intensity / accrual_adj;

        HullBrigoWwrCdsAnalytics {
            clean_risk_free_cds_spread: s_clean,
            systemic_joint_default_intensity: lam_sys,
            wrong_way_cva_spread_haircut: s_clean - s_wwr,
            wrong_way_adjusted_cds_par_spread: s_wwr,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hull_brigo_wwr_cds_invariants_and_ad_parity() {
            let base = HullBrigoWwrCdsSpec {
                reference_entity_hazard_lambda_ref: 0.024_f64,
                protection_seller_hazard_lambda_cp: 0.016,
                default_correlation_rho_wwr: 0.45,
                reference_recovery_rate_r_ref: 0.40,
                counterparty_recovery_rate_r_cp: 0.35,
                cds_spread_volatility_sigma_s: 0.42,
                risk_free_rate_r: 0.035,
                maturity_years_t: 5.0,
            };
            let r0 = evaluate_hull_white_brigo_cds_counterparty_wrong_way_spread(&base);
            assert!(r0.clean_risk_free_cds_spread > r0.wrong_way_adjusted_cds_par_spread);
            assert!(r0.wrong_way_cva_spread_haircut > 0.0);
            assert!(r0.wrong_way_adjusted_cds_par_spread > 0.0);

            let inputs = [
                base.reference_entity_hazard_lambda_ref,
                base.protection_seller_hazard_lambda_cp,
                base.default_correlation_rho_wwr,
                base.reference_recovery_rate_r_ref,
                base.counterparty_recovery_rate_r_cp,
                base.cds_spread_volatility_sigma_s,
                base.risk_free_rate_r,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HullBrigoWwrCdsSpec {
                reference_entity_hazard_lambda_ref: dual_in[0],
                protection_seller_hazard_lambda_cp: dual_in[1],
                default_correlation_rho_wwr: dual_in[2],
                reference_recovery_rate_r_ref: dual_in[3],
                counterparty_recovery_rate_r_cp: dual_in[4],
                cds_spread_volatility_sigma_s: dual_in[5],
                risk_free_rate_r: dual_in[6],
                maturity_years_t: dual_in[7],
            };
            let r_ad = evaluate_hull_white_brigo_cds_counterparty_wrong_way_spread(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| HullBrigoWwrCdsSpec {
                    reference_entity_hazard_lambda_ref: v[0],
                    protection_seller_hazard_lambda_cp: v[1],
                    default_correlation_rho_wwr: v[2],
                    reference_recovery_rate_r_ref: v[3],
                    counterparty_recovery_rate_r_cp: v[4],
                    cds_spread_volatility_sigma_s: v[5],
                    risk_free_rate_r: v[6],
                    maturity_years_t: v[7],
                };
                let rp = evaluate_hull_white_brigo_cds_counterparty_wrong_way_spread(&mk(p));
                let rm = evaluate_hull_white_brigo_cds_counterparty_wrong_way_spread(&mk(m));
                let fd = (rp.wrong_way_adjusted_cds_par_spread
                    - rm.wrong_way_adjusted_cds_par_spread)
                    / (2.0 * eps);
                assert!(
                    r_ad.wrong_way_adjusted_cds_par_spread.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.wrong_way_adjusted_cds_par_spread.d[i] - fd).abs() < 1e-5);
            }
        }
    }
}
pub use hull_brigo_wwr_cds_mod::*;
