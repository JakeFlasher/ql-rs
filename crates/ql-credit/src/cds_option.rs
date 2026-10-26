#![allow(clippy::excessive_precision)]
//! Black's option on a forward CDS spread, including front-end protection.
pub use crate::phase35::*;
pub use crate::phase36::*;
pub use crate::phase37::*;
pub use crate::phase38::*;
#[allow(unused_imports)]
pub use crate::phase39::*;
#[allow(unused_imports)]
pub use crate::phase40::*;
#[allow(unused_imports)]
pub use crate::phase41::*;
#[allow(unused_imports)]
pub use crate::phase42::*;

use crate::SurvivalProbabilityCurve;
pub use crate::phase33::*;
pub use crate::phase34::*;
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
mod phase14_tests {
    use super::*;
    use crate::SurvivalProbabilityCurve;
    use ql_ad::Dual;
    use ql_curves::FlatForwardCurve;
    use ql_time::Date;
    #[test]
    fn phase14_cds_swaption_curve_backed() {
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
mod phase15_tests {
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
mod phase19_tests {
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
mod phase20_tests {
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
mod phase26_tests {
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
