//! Engine 5A: Qian, Hua & Sorensen (2007) & Grinold–Kahn (2000) Multi-Period Information Ratio
//! (`MP-IR`), Alpha Decay Horizon & Turnover-Optimal Trading Speed Engine.
//! Engine 5B: Roncalli (2013) & Meucci (2005) Choueifaty–Coignard Diversification Ratio (`DR`) &
//! Shannon-Entropy Effective Number of Bets (`ENB`) Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct QianMultiPeriodAlphaDecayIrSpec<T: Real = f64> {
    pub initial_ic0: T,
    pub alpha_decay_rate_kappa: T,
    pub breadth_br: T,
    pub proportional_tcost_rate_ctc: T,
    pub active_risk_budget_sigma_a: T,
}

#[derive(Debug, Clone, Copy)]
pub struct QianMultiPeriodAlphaDecayIrReport<T: Real = f64> {
    pub alpha_half_life_years: T,
    pub optimal_rebalancing_speed_eta: T,
    pub realized_multi_period_ic: T,
    pub gross_information_ratio: T,
    pub net_information_ratio: T,
    pub net_active_value_added: T,
}

pub fn evaluate_qian_multi_period_alpha_decay_ir<T: Real>(
    spec: &QianMultiPeriodAlphaDecayIrSpec<T>,
) -> QianMultiPeriodAlphaDecayIrReport<T> {
    let kap = spec.alpha_decay_rate_kappa.max_val(T::from_f64(1e-4));
    let ic0 = spec.initial_ic0.max_val(T::from_f64(1e-4));
    let sqrt_br = spec.breadth_br.max_val(T::one()).sqrt();
    let ctc = spec.proportional_tcost_rate_ctc.max_val(T::from_f64(1e-5));

    let half_life = T::from_f64(std::f64::consts::LN_2) / kap;

    // Optimal trading rate eta^* maximizing IC_0 * sqrt(BR) * eta / (eta + kappa) - c_tc * eta:
    // eta^* = sqrt(kappa * IC_0 * sqrt(BR) / c_tc) - kappa
    let raw_eta = (kap * ic0 * sqrt_br / ctc).sqrt() - kap;
    let eta_star = T::from_f64(0.5) * (raw_eta + (raw_eta * raw_eta + T::from_f64(1e-6)).sqrt());

    let ic_real = ic0 * (eta_star / (eta_star + kap));
    let ir_gross = ic_real * sqrt_br;
    let ir_net = ir_gross - ctc * eta_star;
    let net_va = ir_net * spec.active_risk_budget_sigma_a;

    QianMultiPeriodAlphaDecayIrReport {
        alpha_half_life_years: half_life,
        optimal_rebalancing_speed_eta: eta_star,
        realized_multi_period_ic: ic_real,
        gross_information_ratio: ir_gross,
        net_information_ratio: ir_net,
        net_active_value_added: net_va,
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RoncalliMeucciDiversificationEnbSpec<T: Real = f64> {
    pub weights: [T; 3],
    pub asset_volatilities: [T; 3],
    pub corr_01: T,
    pub corr_02: T,
    pub corr_12: T,
}

#[derive(Debug, Clone, Copy)]
pub struct RoncalliMeucciDiversificationEnbReport<T: Real = f64> {
    pub weighted_average_volatility: T,
    pub portfolio_volatility: T,
    pub choueifaty_diversification_ratio: T,
    pub risk_budget_shares_p: [T; 3],
    pub shannon_entropy_effective_bets_enb: T,
    pub risk_share_partition_error: T,
}

pub fn evaluate_roncalli_meucci_diversification_and_enb<T: Real>(
    spec: &RoncalliMeucciDiversificationEnbSpec<T>,
) -> RoncalliMeucciDiversificationEnbReport<T> {
    let w = spec.weights;
    let sig = spec.asset_volatilities;
    let corr = [
        [T::one(), spec.corr_01, spec.corr_02],
        [spec.corr_01, T::one(), spec.corr_12],
        [spec.corr_02, spec.corr_12, T::one()],
    ];

    let mut weighted_vol = T::zero();
    let mut port_var = T::zero();
    let mut mrc_num = [T::zero(); 3];
    for i in 0..3 {
        weighted_vol += w[i] * sig[i];
        let mut s_i = T::zero();
        for j in 0..3 {
            s_i += sig[i] * sig[j] * corr[i][j] * w[j];
        }
        mrc_num[i] = w[i] * s_i;
        port_var += mrc_num[i];
    }
    let port_vol = port_var.max_val(T::from_f64(1e-12)).sqrt();
    let dr = weighted_vol / port_vol;

    let mut p = [T::zero(); 3];
    let mut p_sum = T::zero();
    let mut neg_ent = T::zero();
    for i in 0..3 {
        let pi = (mrc_num[i] / port_var.max_val(T::from_f64(1e-12))).max_val(T::from_f64(1e-8));
        p[i] = pi;
        p_sum += pi;
        neg_ent -= pi * pi.ln();
    }
    let enb = neg_ent.exp();
    let part_err = (p_sum - T::one()).abs();

    RoncalliMeucciDiversificationEnbReport {
        weighted_average_volatility: weighted_vol,
        portfolio_volatility: port_vol,
        choueifaty_diversification_ratio: dr,
        risk_budget_shares_p: p,
        shannon_entropy_effective_bets_enb: enb,
        risk_share_partition_error: part_err,
    }
}

#[cfg(test)]
mod phase40_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn qian_alpha_decay_and_roncalli_meucci_enb_invariants_and_ad_parity() {
        let s1 = QianMultiPeriodAlphaDecayIrSpec {
            initial_ic0: 0.075,
            alpha_decay_rate_kappa: 2.0,
            breadth_br: 180.0,
            proportional_tcost_rate_ctc: 0.08,
            active_risk_budget_sigma_a: 0.04,
        };
        let r1 = evaluate_qian_multi_period_alpha_decay_ir(&s1);
        assert!(r1.optimal_rebalancing_speed_eta > 0.0);
        assert!(
            r1.net_information_ratio > 0.0 && r1.net_information_ratio < r1.gross_information_ratio
        );

        let h1 = 1e-6;
        let r1_ad = evaluate_qian_multi_period_alpha_decay_ir(&QianMultiPeriodAlphaDecayIrSpec {
            initial_ic0: Dual::<1>::variable(s1.initial_ic0, 0),
            alpha_decay_rate_kappa: Dual::<1>::from_value(s1.alpha_decay_rate_kappa),
            breadth_br: Dual::<1>::from_value(s1.breadth_br),
            proportional_tcost_rate_ctc: Dual::<1>::from_value(s1.proportional_tcost_rate_ctc),
            active_risk_budget_sigma_a: Dual::<1>::from_value(s1.active_risk_budget_sigma_a),
        });
        let mut s1p = s1;
        s1p.initial_ic0 += h1;
        let mut s1m = s1;
        s1m.initial_ic0 -= h1;
        let fd1 = (evaluate_qian_multi_period_alpha_decay_ir(&s1p).net_active_value_added
            - evaluate_qian_multi_period_alpha_decay_ir(&s1m).net_active_value_added)
            / (2.0 * h1);
        assert!((r1_ad.net_active_value_added.d[0] - fd1).abs() < 1e-6);

        let s2 = RoncalliMeucciDiversificationEnbSpec {
            weights: [0.30, 0.35, 0.35],
            asset_volatilities: [0.18, 0.14, 0.12],
            corr_01: 0.25,
            corr_02: 0.10,
            corr_12: 0.30,
        };
        let r2 = evaluate_roncalli_meucci_diversification_and_enb(&s2);
        assert!(r2.choueifaty_diversification_ratio > 1.0);
        assert!(
            r2.shannon_entropy_effective_bets_enb > 1.0
                && r2.shannon_entropy_effective_bets_enb <= 3.0
        );
        assert!(r2.risk_share_partition_error < 1e-12);

        let h = 1e-6;
        let r_ad = evaluate_roncalli_meucci_diversification_and_enb(
            &RoncalliMeucciDiversificationEnbSpec {
                weights: s2.weights.map(Dual::<1>::from_value),
                asset_volatilities: [
                    Dual::<1>::variable(s2.asset_volatilities[0], 0),
                    Dual::<1>::from_value(s2.asset_volatilities[1]),
                    Dual::<1>::from_value(s2.asset_volatilities[2]),
                ],
                corr_01: Dual::<1>::from_value(s2.corr_01),
                corr_02: Dual::<1>::from_value(s2.corr_02),
                corr_12: Dual::<1>::from_value(s2.corr_12),
            },
        );
        let mut sp = s2;
        sp.asset_volatilities[0] += h;
        let mut sm = s2;
        sm.asset_volatilities[0] -= h;
        let fd = (evaluate_roncalli_meucci_diversification_and_enb(&sp)
            .shannon_entropy_effective_bets_enb
            - evaluate_roncalli_meucci_diversification_and_enb(&sm)
                .shannon_entropy_effective_bets_enb)
            / (2.0 * h);
        assert!((r_ad.shannon_entropy_effective_bets_enb.d[0] - fd).abs() < 1e-6);
    }
}
