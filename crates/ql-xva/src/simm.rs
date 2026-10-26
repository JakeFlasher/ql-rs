use ql_ad::Real;

/// The six ISDA SIMM risk classes (in the ordering used by the v2.6
/// cross-risk-class correlation matrix).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
#[repr(usize)]
pub enum SimmRiskClass {
    InterestRate = 0,
    CreditQualifying = 1,
    CreditNonQualifying = 2,
    Equity = 3,
    Commodity = 4,
    Fx = 5,
}

/// Aggregate class margins using the ISDA SIMM v2.6 cross-risk-class matrix.
/// The input ordering is `[IR, CreditQ, CreditNQ, Equity, Commodity, FX]`.
pub fn compute_isda_simm_6_risk_classes<T: Real>(class_margins: &[T; 6]) -> T {
    // ISDA SIMM v2.6, Table 2.2.  The matrix is symmetric with unit diagonal.
    const PSI: [[f64; 6]; 6] = [
        [1.00, 0.04, 0.04, 0.07, 0.37, 0.14],
        [0.04, 1.00, 0.54, 0.70, 0.27, 0.37],
        [0.04, 0.54, 1.00, 0.46, 0.24, 0.15],
        [0.07, 0.70, 0.46, 1.00, 0.35, 0.39],
        [0.37, 0.27, 0.24, 0.35, 1.00, 0.35],
        [0.14, 0.37, 0.15, 0.39, 0.35, 1.00],
    ];
    let mut total = T::zero();
    for i in 0..6 {
        for j in 0..6 {
            total += T::from_f64(PSI[i][j]) * class_margins[i] * class_margins[j];
        }
    }
    // Numerical roundoff can make a nearly-zero portfolio marginally negative
    // for scalar f64.  Keep AD values in-domain without stripping tangents.
    total.max_val(T::zero()).sqrt()
}
#[derive(Clone, Debug)]
pub struct SimmBucket<T: Real> {
    pub sensitivities: Vec<T>,
    pub risk_weights: Vec<f64>,
    pub concentration_threshold: f64,
    pub intra_corr: Vec<Vec<f64>>,
}
pub fn compute_isda_simm<T: Real>(b: &[SimmBucket<T>], corr: &[Vec<f64>]) -> T {
    assert_eq!(
        corr.len(),
        b.len(),
        "inter-bucket correlation dimension mismatch"
    );
    assert!(corr.iter().all(|row| row.len() == b.len()));
    for x in b {
        assert_eq!(x.sensitivities.len(), x.risk_weights.len());
        assert_eq!(x.intra_corr.len(), x.sensitivities.len());
        assert!(
            x.intra_corr
                .iter()
                .all(|row| row.len() == x.sensitivities.len())
        );
        assert!(x.concentration_threshold.is_finite() && x.concentration_threshold > 0.0);
    }
    let mut ks = Vec::new();
    let mut ss = Vec::new();
    for x in b {
        let concentration_ratios: Vec<T> = x
            .sensitivities
            .iter()
            .map(|s| {
                let threshold = x.concentration_threshold.max(1e-12);
                let abs_ratio = s.abs() / threshold;
                // Keep the concentration multiplier in the AD domain.  Using
                // only `s.val()` here would silently drop the derivative of
                // the ISDA concentration adjustment for large sensitivities.
                if abs_ratio.val() > 1.0 {
                    abs_ratio.sqrt()
                } else if abs_ratio.val() == 1.0 {
                    // At the kink use the symmetric finite-difference
                    // subgradient (the mean of the two one-sided slopes).
                    T::one() + (abs_ratio - T::one()) * 0.25
                } else {
                    T::one()
                }
            })
            .collect();
        let ws: Vec<T> = x
            .sensitivities
            .iter()
            .zip(&x.risk_weights)
            .zip(&concentration_ratios)
            .map(|((s, w), concentration)| *s * *concentration * *w)
            .collect();
        let mut k2 = T::zero();
        for i in 0..ws.len() {
            for j in 0..ws.len() {
                // SIMM v2.6 scales off-diagonal correlations by the ratio of
                // concentration factors.  Compute this in T so AD sensitivities
                // remain exact through the concentration adjustment.
                let concentration_corr = if i == j {
                    T::one()
                } else {
                    let lo = concentration_ratios[i].min_val(concentration_ratios[j]);
                    let hi = concentration_ratios[i].max_val(concentration_ratios[j]);
                    if hi.val().abs() > 1e-14 {
                        lo / hi
                    } else {
                        T::one()
                    }
                };
                k2 += T::from_f64(
                    x.intra_corr
                        .get(i)
                        .and_then(|r| r.get(j))
                        .copied()
                        .unwrap_or(if i == j { 1. } else { 0. }),
                ) * concentration_corr
                    * ws[i]
                    * ws[j];
            }
        }
        let k = k2.max_val(T::zero()).sqrt();
        ks.push(k);
        ss.push(
            ws.iter()
                .copied()
                .fold(T::zero(), |a, z| a + z)
                .min_val(k)
                .max_val(-k),
        );
    }
    let mut total = T::zero();
    for i in 0..ks.len() {
        total += ks[i] * ks[i];
        for j in 0..ks.len() {
            if i != j {
                total += T::from_f64(corr.get(i).and_then(|r| r.get(j)).copied().unwrap_or(0.))
                    * ss[i]
                    * ss[j];
            }
        }
    }
    total.max_val(T::zero()).sqrt()
}

#[derive(Clone, Debug)]
pub struct DynamicSimmMvaSpec<T: Real = f64> {
    pub times: Vec<f64>,
    pub bucket_sensitivities_t0: Vec<T>,
    pub bucket_decay_rates: Vec<T>,
    pub risk_weights: Vec<f64>,
    pub concentration_thresholds: Vec<f64>,
    pub inter_bucket_corr: f64,
    pub scheduled_cashflows: Vec<T>,
    pub mpor_years: f64,
    pub margin_funding_spread: T,
    pub hazard_rate: T,
    pub lgd: T,
    pub discount_rate: T,
}

#[derive(Clone, Debug)]
pub struct DynamicSimmMvaReport<T: Real> {
    pub simm_profile: Vec<T>,
    pub mva: T,
    pub base_mpor_cva: T,
    pub cashflow_spike_cva_increment: T,
    pub total_aps_cva: T,
}

pub fn compute_dynamic_simm_mva_and_spike<T: Real>(
    spec: &DynamicSimmMvaSpec<T>,
) -> DynamicSimmMvaReport<T> {
    let n = spec.times.len().min(spec.scheduled_cashflows.len());
    let nb = spec
        .bucket_sensitivities_t0
        .len()
        .min(spec.bucket_decay_rates.len())
        .min(spec.risk_weights.len())
        .min(spec.concentration_thresholds.len());
    if n == 0 || nb == 0 {
        return DynamicSimmMvaReport {
            simm_profile: Vec::new(),
            mva: T::zero(),
            base_mpor_cva: T::zero(),
            cashflow_spike_cva_increment: T::zero(),
            total_aps_cva: T::zero(),
        };
    }
    let t_max = spec.times[n - 1].max(1e-6);
    let mpor = spec.mpor_years.max(1e-4);
    let q99 = 2.3263478740408408; // \Phi^{-1}(0.99)
    let phi0 = 0.3989422804014327; // \phi(0)
    let mut simm_profile = Vec::with_capacity(n);
    let mut mva = T::zero();
    let mut base_cva = T::zero();
    let mut spike_cva = T::zero();
    let mut prev_t = 0.0;

    for k in 0..n {
        let t = spec.times[k].max(prev_t + 1e-6);
        let dt = t - prev_t;
        let rem_frac = (1.0 - t / (t_max + mpor)).max(0.05);
        let mut ws = Vec::with_capacity(nb);
        let mut crs = Vec::with_capacity(nb);
        for b in 0..nb {
            let s_bt = spec.bucket_sensitivities_t0[b]
                * (-spec.bucket_decay_rates[b] * t).exp()
                * rem_frac;
            let thr = spec.concentration_thresholds[b].max(1e-12);
            let ratio = (s_bt * s_bt + T::from_f64(1e-12)).sqrt() / thr;
            let sqrt_ratio = ratio.sqrt();
            let diff_cr = sqrt_ratio - T::one();
            let cr = (T::one() + sqrt_ratio + (diff_cr * diff_cr + T::from_f64(1e-6)).sqrt()) * 0.5;
            crs.push(cr);
            ws.push(s_bt * spec.risk_weights[b] * cr);
        }
        let mut im2 = T::zero();
        for a in 0..nb {
            for b in 0..nb {
                let rho = if a == b {
                    T::one()
                } else {
                    let lo = crs[a].min_val(crs[b]);
                    let hi = crs[a].max_val(crs[b]).max_val(T::from_f64(1e-12));
                    (lo / hi) * spec.inter_bucket_corr
                };
                im2 += ws[a] * ws[b] * rho;
            }
        }
        let im_t = im2.max_val(T::from_f64(1e-16)).sqrt();
        simm_profile.push(im_t);

        let df_surv = (-(spec.discount_rate + spec.hazard_rate) * t).exp();
        mva += spec.margin_funding_spread * im_t * df_surv * dt;

        // Andersen-Pykhtin-Sokol (2017): Brownian bridge MPoR sigma = IM_t / q_99
        let sigma_mpor = im_t / q99;
        let epe_base = sigma_mpor * phi0;
        let cf_k = spec.scheduled_cashflows[k];
        // Unpaid cashflow spike over MPoR window adds Bachelier E[(CF_k + sigma_mpor * Z)^+] - E[(sigma_mpor * Z)^+]
        let z_spike = cf_k / sigma_mpor.max_val(T::from_f64(1e-12));
        let epe_with_spike = cf_k * z_spike.norm_cdf() + sigma_mpor * z_spike.norm_pdf();
        let spike_inc =
            (epe_with_spike - epe_base).max_val(T::zero()) * (mpor / dt.max(mpor)).min(1.0);

        let dp_def = spec.hazard_rate * df_surv * dt;
        base_cva += spec.lgd * epe_base * dp_def;
        spike_cva += spec.lgd * spike_inc * dp_def;
        prev_t = t;
    }

    DynamicSimmMvaReport {
        simm_profile,
        mva,
        base_mpor_cva: base_cva,
        cashflow_spike_cva_increment: spike_cva,
        total_aps_cva: base_cva + spike_cva,
    }
}

#[cfg(test)]
mod domain_suite_31_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn dynamic_simm_mva_and_aps_cashflow_spike_invariants_and_ad() {
        let mut spec = DynamicSimmMvaSpec {
            times: vec![0.5, 1.0, 1.5, 2.0],
            bucket_sensitivities_t0: vec![120.0, 80.0],
            bucket_decay_rates: vec![0.25, 0.30],
            risk_weights: vec![0.50, 0.65],
            concentration_thresholds: vec![100.0, 100.0],
            inter_bucket_corr: 0.35,
            scheduled_cashflows: vec![0.0; 4],
            mpor_years: 10.0 / 252.0,
            margin_funding_spread: 0.015,
            hazard_rate: 0.02,
            lgd: 0.60,
            discount_rate: 0.03,
        };
        let zero_cf = compute_dynamic_simm_mva_and_spike(&spec);
        assert!(zero_cf.mva > 0.0);
        assert!(zero_cf.cashflow_spike_cva_increment.abs() < 1e-12);
        assert!((zero_cf.total_aps_cva - zero_cf.base_mpor_cva).abs() < 1e-12);

        spec.scheduled_cashflows = vec![15.0, 15.0, 15.0, 25.0];
        let with_spike = compute_dynamic_simm_mva_and_spike(&spec);
        assert!(with_spike.cashflow_spike_cva_increment > 0.0);
        assert!(with_spike.total_aps_cva > with_spike.base_mpor_cva);

        let h = 1e-5;
        let eval_s0 = |s0: f64| {
            let mut sp = spec.clone();
            sp.bucket_sensitivities_t0[0] = s0;
            compute_dynamic_simm_mva_and_spike(&sp).mva
        };
        let fd = (eval_s0(120.0 + h) - eval_s0(120.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let spec_d = DynamicSimmMvaSpec {
            times: vec![0.5, 1.0, 1.5, 2.0],
            bucket_sensitivities_t0: vec![D::variable(120.0, 0), D::from_value(80.0)],
            bucket_decay_rates: vec![D::from_value(0.25), D::from_value(0.30)],
            risk_weights: vec![0.50, 0.65],
            concentration_thresholds: vec![100.0, 100.0],
            inter_bucket_corr: 0.35,
            scheduled_cashflows: vec![
                D::from_value(15.0),
                D::from_value(15.0),
                D::from_value(15.0),
                D::from_value(25.0),
            ],
            mpor_years: 10.0 / 252.0,
            margin_funding_spread: D::from_value(0.015),
            hazard_rate: D::from_value(0.02),
            lgd: D::from_value(0.60),
            discount_rate: D::from_value(0.03),
        };
        let ad = compute_dynamic_simm_mva_and_spike(&spec_d).mva;
        assert!(ad.d[0] > 0.0);
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

mod evaluate_gregory_bilateral_simm_vs_ccp_mva_module {
    #![allow(
        clippy::needless_range_loop,
        clippy::assign_op_pattern,
        clippy::approx_constant
    )]
    use ql_ad::Real;

    #[derive(Clone, Copy, Debug)]
    pub struct GregoryBilateralSimmVsCcpMvaSpec<T: Real = f64> {
        pub portfolio_annual_volatility: T,
        pub bilateral_mpor_years: T,
        pub ccp_mpor_years: T,
        pub exposure_decay_kappa: T,
        pub im_funding_spread: T,
        pub cpty_hazard_rate: T,
        pub cpty_lgd: T,
        pub risk_free_rate: T,
        pub ccp_default_fund_contribution: T,
        pub capital_hurdle_rate: T,
        pub maturity_years: T,
    }

    #[derive(Clone, Copy, Debug)]
    pub struct GregoryBilateralSimmVsCcpMvaReport<T: Real = f64> {
        pub bilateral_simm_initial_margin: T,
        pub ccp_initial_margin: T,
        pub mpor_sqrt_scaling_ratio: T,
        pub bilateral_mva_charge: T,
        pub bilateral_residual_mpor_cva: T,
        pub total_bilateral_uncleared_xva: T,
        pub ccp_mva_charge: T,
        pub ccp_default_fund_kva: T,
        pub total_ccp_cleared_xva: T,
        pub central_clearing_xva_savings: T,
    }

    pub fn evaluate_gregory_bilateral_simm_vs_ccp_mva<T: Real>(
        spec: &GregoryBilateralSimmVsCcpMvaSpec<T>,
    ) -> GregoryBilateralSimmVsCcpMvaReport<T> {
        let z99 = T::from_f64(2.326_347_874_040_840_8);
        let tau_b = spec.bilateral_mpor_years.max_val(T::from_f64(1e-6));
        let tau_c = spec.ccp_mpor_years.max_val(T::from_f64(1e-6));
        let im_b = z99 * spec.portfolio_annual_volatility * tau_b.sqrt();
        let im_c = z99 * spec.portfolio_annual_volatility * tau_c.sqrt();
        let mpor_ratio = (tau_b / tau_c).sqrt();

        let tot_disc = (spec.risk_free_rate + spec.cpty_hazard_rate + spec.exposure_decay_kappa)
            .max_val(T::from_f64(1e-6));
        let mat = spec.maturity_years.max_val(T::from_f64(1e-6));
        let ann = (T::one() - (-tot_disc * mat).exp()) / tot_disc;

        let mva_b = spec.im_funding_spread * im_b * ann;
        let inv_sqrt_2pi = T::from_f64(0.398_942_280_401_432_7);
        let res_cva_b = spec.cpty_lgd
            * spec.cpty_hazard_rate
            * (spec.portfolio_annual_volatility * tau_b.sqrt() * inv_sqrt_2pi)
            * ann;
        let tot_b = mva_b + res_cva_b;

        let mva_c = spec.im_funding_spread * im_c * ann;
        let df_kva_c = spec.capital_hurdle_rate * spec.ccp_default_fund_contribution * ann;
        let tot_c = mva_c + df_kva_c;

        GregoryBilateralSimmVsCcpMvaReport {
            bilateral_simm_initial_margin: im_b,
            ccp_initial_margin: im_c,
            mpor_sqrt_scaling_ratio: mpor_ratio,
            bilateral_mva_charge: mva_b,
            bilateral_residual_mpor_cva: res_cva_b,
            total_bilateral_uncleared_xva: tot_b,
            ccp_mva_charge: mva_c,
            ccp_default_fund_kva: df_kva_c,
            total_ccp_cleared_xva: tot_c,
            central_clearing_xva_savings: tot_b - tot_c,
        }
    }

    #[cfg(test)]
    mod domain_suite_49_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gregory_bilateral_simm_vs_ccp_mva_invariants_and_ad_parity() {
            let spec = GregoryBilateralSimmVsCcpMvaSpec {
                portfolio_annual_volatility: 15.0,
                bilateral_mpor_years: 10.0 / 250.0,
                ccp_mpor_years: 5.0 / 250.0,
                exposure_decay_kappa: 0.12,
                im_funding_spread: 0.0085,
                cpty_hazard_rate: 0.015,
                cpty_lgd: 0.60,
                risk_free_rate: 0.035,
                ccp_default_fund_contribution: 0.35,
                capital_hurdle_rate: 0.10,
                maturity_years: 5.0,
            };
            let rep = evaluate_gregory_bilateral_simm_vs_ccp_mva(&spec);
            assert!(rep.bilateral_simm_initial_margin > rep.ccp_initial_margin);
            assert!(rep.bilateral_mva_charge > rep.ccp_mva_charge);

            let h = 1e-6;
            let r_ad =
                evaluate_gregory_bilateral_simm_vs_ccp_mva(&GregoryBilateralSimmVsCcpMvaSpec {
                    portfolio_annual_volatility: Dual::<1>::variable(
                        spec.portfolio_annual_volatility,
                        0,
                    ),
                    bilateral_mpor_years: Dual::<1>::from_value(spec.bilateral_mpor_years),
                    ccp_mpor_years: Dual::<1>::from_value(spec.ccp_mpor_years),
                    exposure_decay_kappa: Dual::<1>::from_value(spec.exposure_decay_kappa),
                    im_funding_spread: Dual::<1>::from_value(spec.im_funding_spread),
                    cpty_hazard_rate: Dual::<1>::from_value(spec.cpty_hazard_rate),
                    cpty_lgd: Dual::<1>::from_value(spec.cpty_lgd),
                    risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                    ccp_default_fund_contribution: Dual::<1>::from_value(
                        spec.ccp_default_fund_contribution,
                    ),
                    capital_hurdle_rate: Dual::<1>::from_value(spec.capital_hurdle_rate),
                    maturity_years: Dual::<1>::from_value(spec.maturity_years),
                });
            let mut sp = spec;
            sp.portfolio_annual_volatility += h;
            let mut sm = spec;
            sm.portfolio_annual_volatility -= h;
            let fd = (evaluate_gregory_bilateral_simm_vs_ccp_mva(&sp).central_clearing_xva_savings
                - evaluate_gregory_bilateral_simm_vs_ccp_mva(&sm).central_clearing_xva_savings)
                / (2.0 * h);
            assert!((r_ad.central_clearing_xva_savings.d[0] - fd).abs() < 1e-5);
        }
    }
}
pub use evaluate_gregory_bilateral_simm_vs_ccp_mva_module::*;
