#![allow(clippy::needless_range_loop)]

pub use crate::phase30::*;
pub use crate::phase31::{
    BrockHommesAdaptiveBeliefReport, BrockHommesAdaptiveBeliefSpec, DgtwTreynorMazuyTimingReport,
    DgtwTreynorMazuyTimingSpec, evaluate_brock_hommes_adaptive_belief_switching,
    evaluate_dgtw_characteristic_and_treynor_mazuy_timing,
};
pub use crate::phase32::*;
use ql_ad::Real;

#[derive(Clone, Debug)]
pub struct BussetiBoydRiskConstrainedKellySpec<T: Real = f64> {
    pub excess_return_means: Vec<T>,
    pub covariance: Vec<Vec<T>>,
    pub risk_free_rate: T,
    pub drawdown_floor_alpha: T,
    pub max_drawdown_prob_beta: T,
    pub specialist_immediate_capital: T,
    pub slow_moving_external_capital: T,
    pub capital_arrival_rate: T,
    pub dislocation_horizon_years: T,
    pub forced_liquidation_notional: T,
}
#[derive(Clone, Debug)]
pub struct BussetiBoydRiskConstrainedKellyReport<T: Real = f64> {
    pub full_kelly_weights: Vec<T>,
    pub risk_constrained_kelly_weights: Vec<T>,
    pub drawdown_risk_aversion_lambda: T,
    pub kelly_fraction_kappa: T,
    pub full_kelly_expected_log_growth: T,
    pub constrained_kelly_expected_log_growth: T,
    pub constrained_portfolio_volatility: T,
    pub initial_fire_sale_dislocation_bps: T,
    pub horizon_fire_sale_discount_bps: T,
    pub half_life_capital_recovery_days: T,
}
pub fn solve_busseti_boyd_risk_constrained_kelly<T: Real>(
    s: &BussetiBoydRiskConstrainedKellySpec<T>,
) -> BussetiBoydRiskConstrainedKellyReport<T> {
    let n = s.excess_return_means.len();
    let w = solve_dense_alloc(&s.covariance, &s.excess_return_means);
    let lam = (T::one() / s.max_drawdown_prob_beta.max_val(T::from_f64(1e-12))).ln()
        / (T::one() / s.drawdown_floor_alpha.max_val(T::from_f64(1e-12))).ln();
    let kap = (T::from_f64(2.0) / (T::one() + lam))
        .min_val(T::one())
        .max_val(T::zero());
    let wc = w.iter().map(|&x| x * kap).collect::<Vec<_>>();
    let quad = |x: &[T]| {
        (0..n).fold(T::zero(), |acc, i| {
            acc + (0..n).fold(T::zero(), |a, j| {
                a + x[i]
                    * s.covariance
                        .get(i)
                        .and_then(|r| r.get(j))
                        .copied()
                        .unwrap_or(T::zero())
                    * x[j]
            })
        })
    };
    let dot = |x: &[T]| {
        x.iter()
            .zip(s.excess_return_means.iter())
            .map(|(&a, &b)| a * b)
            .fold(T::zero(), |a, b| a + b)
    };
    let vol_full = quad(&w).max_val(T::from_f64(1e-16)).sqrt();
    let mu_full = dot(&w);
    let gf = s.risk_free_rate + mu_full - T::from_f64(0.5) * vol_full * vol_full;
    let gc = s.risk_free_rate + dot(&wc) - T::from_f64(0.5) * quad(&wc);
    let vol = quad(&wc).max_val(T::from_f64(1e-16)).sqrt();
    let eff_cap = s.specialist_immediate_capital.max_val(T::from_f64(1e-8))
        + s.slow_moving_external_capital.max_val(T::zero())
            * (T::one()
                - (-s.capital_arrival_rate.abs() * s.dislocation_horizon_years.max_val(T::zero()))
                    .exp());
    let dis = (s.forced_liquidation_notional * vol.powf(2.0) / eff_cap * T::from_f64(1e4)).abs();
    BussetiBoydRiskConstrainedKellyReport {
        full_kelly_weights: w,
        risk_constrained_kelly_weights: wc,
        drawdown_risk_aversion_lambda: lam,
        kelly_fraction_kappa: kap,
        full_kelly_expected_log_growth: gf,
        constrained_kelly_expected_log_growth: gc,
        constrained_portfolio_volatility: vol,
        initial_fire_sale_dislocation_bps: dis,
        horizon_fire_sale_discount_bps: dis
            * (T::one()
                - (-s.capital_arrival_rate.abs() * s.dislocation_horizon_years.max_val(T::zero()))
                    .exp())
            .max_val(T::zero()),
        half_life_capital_recovery_days: T::from_f64(2.0_f64.ln())
            / s.capital_arrival_rate.max_val(T::from_f64(1e-8)),
    }
}

use thiserror::Error;
fn solve_dense_alloc<T: Real>(a: &[Vec<T>], b: &[T]) -> Vec<T> {
    let n = b.len();
    let mut m = vec![vec![T::zero(); n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = a[i][j];
        }
        m[i][n] = b[i];
    }
    for i in 0..n {
        let mut p = i;
        for r in i + 1..n {
            if m[r][i].abs().val() > m[p][i].abs().val() {
                p = r;
            }
        }
        m.swap(i, p);
        let d = if m[i][i].abs().val() < 1e-14 {
            if m[i][i].val() < 0.0 {
                T::from_f64(-1e-14)
            } else {
                T::from_f64(1e-14)
            }
        } else {
            m[i][i]
        };
        for j in i..=n {
            m[i][j] /= d;
        }
        let piv = m[i].clone();
        for r in 0..n {
            if r != i {
                let f = m[r][i];
                for j in i..=n {
                    m[r][j] -= f * piv[j];
                }
            }
        }
    }
    (0..n).map(|i| m[i][n]).collect()
}
#[derive(Clone, Debug)]
pub struct GrinoldKahnClarkeFundamentalLawSpec<T: Real = f64> {
    pub standardized_scores_z: Vec<T>,
    pub residual_volatilities: Vec<T>,
    pub benchmark_weights: Vec<T>,
    pub information_coefficient_ic: T,
    pub rebalances_per_year: usize,
    pub signal_autocorrelation_rho: T,
    pub target_tracking_error: T,
    pub max_active_weight_bound: T,
}
#[derive(Clone, Debug)]
pub struct GrinoldKahnClarkeFundamentalLawReport<T: Real = f64> {
    pub nominal_breadth_br: T,
    pub autocorrelation_adjusted_breadth_br_eff: T,
    pub refined_alphas: Vec<T>,
    pub unconstrained_active_weights: Vec<T>,
    pub constrained_active_weights: Vec<T>,
    pub transfer_coefficient_tc: T,
    pub unconstrained_information_ratio: T,
    pub constrained_information_ratio: T,
    pub expected_active_return: T,
    pub optimal_active_risk_aversion_lambda: T,
}
pub fn evaluate_grinold_kahn_clarke_fundamental_law<T: Real>(
    s: &GrinoldKahnClarkeFundamentalLawSpec<T>,
) -> GrinoldKahnClarkeFundamentalLawReport<T> {
    let n = s
        .standardized_scores_z
        .len()
        .min(s.residual_volatilities.len())
        .min(s.benchmark_weights.len());
    let rebs = T::from_f64(s.rebalances_per_year.max(1) as f64);
    let br = T::from_f64(n as f64) * rebs;
    let rho = s.signal_autocorrelation_rho;
    let bre = br * (T::one() - rho) / (T::one() + rho).max_val(T::from_f64(1e-8));
    let mut alpha = Vec::with_capacity(n);
    let mut wu = Vec::with_capacity(n);
    let mut var_raw = T::zero();
    for i in 0..n {
        let sig = s.residual_volatilities[i].abs().max_val(T::from_f64(1e-8));
        let a = s.information_coefficient_ic * sig * s.standardized_scores_z[i];
        alpha.push(a);
        let w = a / (sig * sig);
        wu.push(w);
        var_raw += w * w * sig * sig;
    }
    let te_target = s.target_tracking_error.max_val(T::from_f64(1e-8));
    let norm_u = var_raw.sqrt().max_val(T::from_f64(1e-8));
    let scale = te_target / norm_u;
    for x in &mut wu {
        *x *= scale;
    }
    let bound = s.max_active_weight_bound.abs().max_val(T::from_f64(1e-8));
    let smooth_pos = |v: T| T::from_f64(0.5) * (v + (v * v + T::from_f64(1e-24)).sqrt());
    let mut wc = Vec::with_capacity(n);
    let mut vc = T::zero();
    for i in 0..n {
        let lo = -s.benchmark_weights[i].abs().min_val(bound);
        let hi = bound;
        let clamped = lo + smooth_pos(wu[i] - lo) - smooth_pos(wu[i] - hi);
        let sig = s.residual_volatilities[i].abs().max_val(T::from_f64(1e-8));
        wc.push(clamped);
        vc += clamped * clamped * sig * sig;
    }
    let scale_c = te_target / vc.sqrt().max_val(T::from_f64(1e-8));
    let mut num_tc = T::zero();
    let mut den_a = T::zero();
    let mut den_w = T::zero();
    let mut rc = T::zero();
    for i in 0..n {
        wc[i] *= scale_c;
        let sig = s.residual_volatilities[i].abs().max_val(T::from_f64(1e-8));
        let ai_std = alpha[i] / sig;
        let wi_std = wc[i] * sig;
        num_tc += ai_std * wi_std;
        den_a += ai_std * ai_std;
        den_w += wi_std * wi_std;
        rc += alpha[i] * wc[i];
    }
    let tc = (num_tc / (den_a.sqrt() * den_w.sqrt()).max_val(T::from_f64(1e-12)))
        .max_val(T::zero())
        .min_val(T::one());
    let iru = s.information_coefficient_ic * bre.sqrt();
    let irc = tc * iru;
    let lam = norm_u / (T::from_f64(2.0) * te_target);
    GrinoldKahnClarkeFundamentalLawReport {
        nominal_breadth_br: br,
        autocorrelation_adjusted_breadth_br_eff: bre,
        refined_alphas: alpha,
        unconstrained_active_weights: wu,
        constrained_active_weights: wc,
        transfer_coefficient_tc: tc,
        unconstrained_information_ratio: iru,
        constrained_information_ratio: irc,
        expected_active_return: rc * rebs.sqrt(),
        optimal_active_risk_aversion_lambda: lam,
    }
}

#[derive(Clone, Debug)]
pub struct LopezDePradoNcoDenoisedSpec<T: Real = f64> {
    pub covariance: Vec<Vec<T>>,
    pub expected_returns: Vec<T>,
    pub sample_ratio_q: T,
    pub denoising_alpha: T,
    pub num_clusters: usize,
}
#[derive(Clone, Debug)]
pub struct LopezDePradoNcoDenoisedReport<T: Real = f64> {
    pub marchenko_pastur_lambda_max: T,
    pub num_signal_eigenvalues: usize,
    pub raw_condition_number: T,
    pub denoised_condition_number: T,
    pub denoised_covariance: Vec<Vec<T>>,
    pub detoned_correlation_matrix: Vec<Vec<T>>,
    pub intra_cluster_weights: Vec<T>,
    pub inter_cluster_allocations: Vec<T>,
    pub nco_optimal_weights: Vec<T>,
    pub standard_markowitz_weights: Vec<T>,
    pub nco_portfolio_volatility: T,
}
pub fn solve_lopez_de_prado_denoised_nco_allocation<T: Real>(
    s: &LopezDePradoNcoDenoisedSpec<T>,
) -> LopezDePradoNcoDenoisedReport<T> {
    let n = s.covariance.len();
    assert!(n > 0);
    let mut vol = vec![T::zero(); n];
    let mut corr = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        vol[i] = s.covariance[i][i].max_val(T::from_f64(1e-12)).sqrt();
    }
    for i in 0..n {
        for j in 0..n {
            corr[i][j] = s.covariance[i][j] / (vol[i] * vol[j]);
        }
    }
    let (eig, vecs) = jacobi_eigen_spd(&corr);
    let q = s.sample_ratio_q.max_val(T::from_f64(1.01));
    let lmax = (T::one() + (T::one() / q).sqrt()).powf(2.0);
    let sig = eig.iter().filter(|x| x.val() > lmax.val()).count();
    let noise = eig
        .iter()
        .filter(|x| x.val() <= lmax.val())
        .copied()
        .fold(T::zero(), |a, v| a + v)
        / T::from_f64((n - sig).max(1) as f64);
    let alpha = s.denoising_alpha.max_val(T::zero()).min_val(T::one());
    let mut deig = eig.clone();
    for i in sig..n {
        deig[i] = alpha * eig[i] + (T::one() - alpha) * noise;
    }
    let mut dc = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                dc[i][j] += vecs[i][k] * deig[k] * vecs[j][k];
            }
        }
    }
    let dc_diag: Vec<T> = (0..n)
        .map(|i| dc[i][i].max_val(T::from_f64(1e-12)).sqrt())
        .collect();
    for i in 0..n {
        for j in 0..n {
            dc[i][j] /= dc_diag[i] * dc_diag[j];
        }
    }
    let mut ds = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            ds[i][j] = vol[i] * vol[j] * dc[i][j];
        }
    }
    let rawc = eig[0] / eig[n - 1].max_val(T::from_f64(1e-12));
    let cond = deig[0] / deig[n - 1].max_val(T::from_f64(1e-12));
    let mut det = dc.clone();
    for i in 0..n {
        for j in 0..n {
            det[i][j] -= deig[0] * vecs[i][0] * vecs[j][0];
        }
    }
    let det_diag: Vec<T> = (0..n)
        .map(|i| det[i][i].abs().max_val(T::from_f64(1e-12)).sqrt())
        .collect();
    for i in 0..n {
        for j in 0..n {
            det[i][j] /= det_diag[i] * det_diag[j];
        }
    }
    let mut mark_a = s.covariance.clone();
    for i in 0..n {
        mark_a[i][i] += T::from_f64(1e-8);
    }
    let zraw = solve_dense_alloc(&mark_a, &s.expected_returns);
    let zsum = zraw.iter().copied().fold(T::zero(), |a, v| a + v);
    let zden = if zsum.abs().val() < 1e-12 {
        T::from_f64(1e-12)
    } else {
        zsum
    };
    let mw: Vec<T> = zraw.iter().map(|v| *v / zden).collect();
    let k = s.num_clusters.clamp(1, n);
    let mut intra = vec![T::zero(); n];
    let ranges: Vec<(usize, usize)> = (0..k).map(|c| (c * n / k, (c + 1) * n / k)).collect();
    for &(lo, hi) in &ranges {
        let sz = hi - lo;
        let mut sub = vec![vec![T::zero(); sz]; sz];
        let mut mu = vec![T::zero(); sz];
        for i in 0..sz {
            mu[i] = s.expected_returns[lo + i];
            for j in 0..sz {
                sub[i][j] = ds[lo + i][lo + j];
            }
            sub[i][i] += T::from_f64(1e-8);
        }
        let u = solve_dense_alloc(&sub, &mu);
        let us = u.iter().copied().fold(T::zero(), |a, v| a + v);
        let ud = if us.abs().val() < 1e-12 {
            T::from_f64(1e-12)
        } else {
            us
        };
        for i in 0..sz {
            intra[lo + i] = u[i] / ud;
        }
    }
    let mut omega = vec![vec![T::zero(); k]; k];
    for aidx in 0..k {
        let (alo, ahi) = ranges[aidx];
        for bidx in 0..k {
            let (blo, bhi) = ranges[bidx];
            for i in alo..ahi {
                for j in blo..bhi {
                    omega[aidx][bidx] += intra[i] * ds[i][j] * intra[j];
                }
            }
        }
    }
    for i in 0..k {
        omega[i][i] += T::from_f64(1e-8);
    }
    let c_rhs = vec![T::one(); k];
    let eraw = solve_dense_alloc(&omega, &c_rhs);
    let esum = eraw.iter().copied().fold(T::zero(), |a, v| a + v);
    let eden = if esum.abs().val() < 1e-12 {
        T::from_f64(1e-12)
    } else {
        esum
    };
    let eta: Vec<T> = eraw.iter().map(|v| *v / eden).collect();
    let mut w = vec![T::zero(); n];
    for i in 0..n {
        let c = ranges
            .iter()
            .position(|&(lo, hi)| i >= lo && i < hi)
            .unwrap_or(0);
        w[i] = eta[c] * intra[i];
    }
    let mut pv = T::zero();
    for i in 0..n {
        for j in 0..n {
            pv += w[i] * ds[i][j] * w[j];
        }
    }
    LopezDePradoNcoDenoisedReport {
        marchenko_pastur_lambda_max: lmax,
        num_signal_eigenvalues: sig,
        raw_condition_number: rawc,
        denoised_condition_number: cond.min_val(rawc),
        denoised_covariance: ds,
        detoned_correlation_matrix: det,
        intra_cluster_weights: intra,
        inter_cluster_allocations: eta,
        nco_optimal_weights: w,
        standard_markowitz_weights: mw,
        nco_portfolio_volatility: pv.max_val(T::zero()).sqrt(),
    }
}

#[derive(Clone, Debug)]
pub struct MeucciIdzorekBlackLittermanSpec<T: Real = f64> {
    pub covariance: Vec<Vec<T>>,
    pub market_weights: Vec<T>,
    pub risk_aversion_lambda: T,
    pub tau: T,
    pub view_matrix_p: Vec<Vec<T>>,
    pub view_returns_q: Vec<T>,
    pub view_confidences: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct MeucciIdzorekBlackLittermanReport<T: Real = f64> {
    pub implied_equilibrium_returns: Vec<T>,
    pub calibrated_omega_diagonal: Vec<T>,
    pub posterior_expected_returns: Vec<T>,
    pub posterior_return_covariance: Vec<Vec<T>>,
    pub optimal_posterior_weights: Vec<T>,
    pub active_weight_tilts: Vec<T>,
    pub mahalanobis_plausibility_squared: T,
    pub view_weight_sensitivities_l2: Vec<T>,
    pub grinold_kahn_transfer_coefficient: T,
    pub meucci_ellipsoidal_shrinkage_factor: T,
}
fn invm<T: Real>(a: &[Vec<T>]) -> Vec<Vec<T>> {
    let n = a.len();
    let mut x = vec![vec![T::zero(); 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            x[i][j] = a[i][j];
        }
        x[i][n + i] = T::one();
    }
    for i in 0..n {
        let mut p = i;
        for r in i + 1..n {
            if x[r][i].abs().val() > x[p][i].abs().val() {
                p = r;
            }
        }
        x.swap(i, p);
        let d = if x[i][i].val() >= 0.0 {
            x[i][i].max_val(T::from_f64(1e-14))
        } else {
            x[i][i].min_val(T::from_f64(-1e-14))
        };
        for j in 0..2 * n {
            x[i][j] /= d;
        }
        for r in 0..n {
            if r != i {
                let f = x[r][i];
                let piv = x[i].clone();
                for j in 0..2 * n {
                    x[r][j] -= f * piv[j];
                }
            }
        }
    }
    x.into_iter().map(|r| r[n..].to_vec()).collect()
}
pub fn solve_meucci_idzorek_black_litterman_and_plausibility<T: Real>(
    s: &MeucciIdzorekBlackLittermanSpec<T>,
) -> MeucciIdzorekBlackLittermanReport<T> {
    let n = s.market_weights.len();
    let k = s.view_returns_q.len();
    let mut pi = vec![T::zero(); n];
    for i in 0..n {
        for j in 0..n {
            pi[i] += s.risk_aversion_lambda * s.covariance[i][j] * s.market_weights[j];
        }
    }
    let mut ts = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            ts[i][j] = s.tau * s.covariance[i][j];
        }
    }
    let mut om = vec![T::zero(); k];
    for a in 0..k {
        let mut v = T::zero();
        for i in 0..n {
            for j in 0..n {
                v += s.view_matrix_p[a][i] * ts[i][j] * s.view_matrix_p[a][j];
            }
        }
        let c = s
            .view_confidences
            .get(a)
            .copied()
            .unwrap_or(T::from_f64(0.5))
            .max_val(T::from_f64(1e-6))
            .min_val(T::from_f64(1.0 - 1e-6));
        om[a] = (T::one() - c) / c * v;
    }
    let mut mtx = vec![vec![T::zero(); k]; k];
    for a in 0..k {
        for b in 0..k {
            for i in 0..n {
                for j in 0..n {
                    mtx[a][b] += s.view_matrix_p[a][i] * ts[i][j] * s.view_matrix_p[b][j];
                }
            }
            if a == b {
                mtx[a][b] += om[a];
            }
        }
    }
    let mi = invm(&mtx);
    let mut kap = vec![vec![T::zero(); k]; n];
    for i in 0..n {
        for b in 0..k {
            for a in 0..k {
                for j in 0..n {
                    kap[i][b] += ts[i][j] * s.view_matrix_p[a][j] * mi[a][b];
                }
            }
        }
    }
    let mut mu = pi.clone();
    for i in 0..n {
        for b in 0..k {
            let mut q = s.view_returns_q[b];
            for j in 0..n {
                q -= s.view_matrix_p[b][j] * pi[j];
            }
            mu[i] += kap[i][b] * q;
        }
    }
    let mut post = s.covariance.clone();
    for i in 0..n {
        for j in 0..n {
            post[i][j] *= T::one() + s.tau;
            for a in 0..k {
                for b in 0..n {
                    post[i][j] -= kap[i][a] * s.view_matrix_p[a][b] * ts[b][j];
                }
            }
        }
    }
    let invs = invm(&s.covariance);
    let mut w = vec![T::zero(); n];
    for i in 0..n {
        for j in 0..n {
            w[i] += invs[i][j] * mu[j] / s.risk_aversion_lambda.max_val(T::from_f64(1e-12));
        }
    }
    let active = w
        .iter()
        .zip(s.market_weights.iter())
        .map(|(a, b)| *a - *b)
        .collect::<Vec<_>>();
    let mut alpha = vec![T::zero(); n];
    for i in 0..n {
        alpha[i] = mu[i] - pi[i];
    }
    let mut maha = T::zero();
    for i in 0..n {
        for j in 0..n {
            maha += alpha[i] / s.tau.max_val(T::from_f64(1e-12)) * invs[i][j] * alpha[j];
        }
    }
    let sens = (0..k)
        .map(|a| {
            let mut z = T::zero();
            for i in 0..n {
                let mut dw = T::zero();
                for j in 0..n {
                    dw += invs[i][j] * kap[j][a];
                }
                z += (dw / s.risk_aversion_lambda.max_val(T::from_f64(1e-12))).powf(2.0);
            }
            z.sqrt()
        })
        .collect();
    let mut proj = w.clone();
    for x in &mut proj {
        *x = x.max_val(T::zero());
    }
    let sumw: T = proj.iter().copied().fold(T::zero(), |a, b| a + b);
    let summ: T = s
        .market_weights
        .iter()
        .copied()
        .fold(T::zero(), |a, b| a + b);
    for x in &mut proj {
        *x *= summ / sumw.max_val(T::from_f64(1e-12));
    }
    let mut num = T::zero();
    let mut den1 = T::zero();
    let mut den2 = T::zero();
    for i in 0..n {
        num += (proj[i] - s.market_weights[i]) * alpha[i];
        for j in 0..n {
            den1 += (proj[i] - s.market_weights[i])
                * s.covariance[i][j]
                * (proj[j] - s.market_weights[j]);
            den2 += alpha[i] * invs[i][j] * alpha[j];
        }
    }
    let tc = (num / (den1.max_val(T::from_f64(1e-16)) * den2.max_val(T::from_f64(1e-16))).sqrt())
        .max_val(T::from_f64(-1.0))
        .min_val(T::one());
    let shrink = (T::one() - s.tau.sqrt() / (den2 + T::from_f64(1e-6)).sqrt()).max_val(T::zero());
    MeucciIdzorekBlackLittermanReport {
        implied_equilibrium_returns: pi,
        calibrated_omega_diagonal: om,
        posterior_expected_returns: mu,
        posterior_return_covariance: post,
        optimal_posterior_weights: w,
        active_weight_tilts: active,
        mahalanobis_plausibility_squared: maha,
        view_weight_sensitivities_l2: sens,
        grinold_kahn_transfer_coefficient: tc,
        meucci_ellipsoidal_shrinkage_factor: shrink,
    }
}

#[derive(Clone, Debug)]
pub struct PaleologoFactorNeutralAlphaSpec<T: Real = f64> {
    pub raw_alphas: Vec<T>,
    pub factor_loadings_b: Vec<Vec<T>>,
    pub factor_covariance_omega: Vec<Vec<T>>,
    pub sample_idio_variances: Vec<T>,
    pub idio_shrinkage_xi: T,
    pub idio_risk_aversion: T,
    pub factor_risk_aversion: T,
    pub factor_crowding_scores: Vec<T>,
    pub crowding_penalty_lambda: T,
}
#[derive(Clone, Debug)]
pub struct PaleologoFactorNeutralAlphaReport<T: Real = f64> {
    pub orthogonal_alphas: Vec<T>,
    pub shrunk_idio_variances: Vec<T>,
    pub pure_factor_neutral_weights: Vec<T>,
    pub crowding_penalized_weights: Vec<T>,
    pub residual_factor_exposures: Vec<T>,
    pub idiosyncratic_information_ratio: T,
    pub proportion_idiosyncratic_variance_pvi: T,
    pub expected_active_return: T,
}
fn paleo_solve<T: Real>(mut a: Vec<Vec<T>>, mut b: Vec<T>) -> Vec<T> {
    let n = b.len();
    for i in 0..n {
        let p = (i..n)
            .max_by(|&x, &y| {
                a[x][i]
                    .val()
                    .abs()
                    .partial_cmp(&a[y][i].val().abs())
                    .unwrap()
            })
            .unwrap();
        a.swap(i, p);
        b.swap(i, p);
        let d = if a[i][i].val() >= 0.0 {
            a[i][i].max_val(T::from_f64(1e-12))
        } else {
            a[i][i].min_val(T::from_f64(-1e-12))
        };
        let row = a[i].clone();
        let rhs = b[i];
        for j in i + 1..n {
            let z = a[j][i] / d;
            for k in i..n {
                a[j][k] -= z * row[k];
            }
            b[j] -= z * rhs;
        }
    }
    let mut x = vec![T::zero(); n];
    for i in (0..n).rev() {
        let mut z = b[i];
        for j in i + 1..n {
            z -= a[i][j] * x[j];
        }
        let d = if a[i][i].val() >= 0.0 {
            a[i][i].max_val(T::from_f64(1e-12))
        } else {
            a[i][i].min_val(T::from_f64(-1e-12))
        };
        x[i] = z / d;
    }
    x
}

pub fn solve_paleologo_factor_neutral_portfolio<T: Real>(
    s: &PaleologoFactorNeutralAlphaSpec<T>,
) -> PaleologoFactorNeutralAlphaReport<T> {
    let n = s.raw_alphas.len();
    let k = s.factor_covariance_omega.len();
    let xi = s.idio_shrinkage_xi.max_val(T::zero()).min_val(T::one());
    let mean = s
        .sample_idio_variances
        .iter()
        .copied()
        .fold(T::zero(), |a, b| a + b)
        / n.max(1) as f64;
    let d: Vec<T> = s
        .sample_idio_variances
        .iter()
        .map(|x| *x * (T::one() - xi) + xi * mean)
        .collect();
    let mut g = vec![vec![T::zero(); k]; k];
    let mut ba = vec![T::zero(); k];
    for i in 0..n {
        let di = T::one() / d[i].max_val(T::from_f64(1e-12));
        for p in 0..k {
            let bp = s
                .factor_loadings_b
                .get(i)
                .and_then(|r| r.get(p))
                .copied()
                .unwrap_or(T::zero());
            ba[p] += bp * di * s.raw_alphas[i];
            for q in 0..k {
                let bq = s
                    .factor_loadings_b
                    .get(i)
                    .and_then(|r| r.get(q))
                    .copied()
                    .unwrap_or(T::zero());
                g[p][q] += bp * di * bq;
            }
        }
    }
    let beta_alpha = paleo_solve(g.clone(), ba.clone());
    let orthogonal_alphas: Vec<T> = (0..n)
        .map(|i| {
            s.raw_alphas[i]
                - (0..k)
                    .map(|p| {
                        s.factor_loadings_b
                            .get(i)
                            .and_then(|r| r.get(p))
                            .copied()
                            .unwrap_or(T::zero())
                            * beta_alpha[p]
                    })
                    .fold(T::zero(), |a, b| a + b)
        })
        .collect();
    let gamma_i = s.idio_risk_aversion.max_val(T::from_f64(1e-8));
    let pure_factor_neutral_weights: Vec<T> = orthogonal_alphas
        .iter()
        .enumerate()
        .map(|(i, x)| *x / d[i].max_val(T::from_f64(1e-12)) / gamma_i)
        .collect();
    let gamma_ratio = s.factor_risk_aversion / gamma_i;
    let mut omega_tilde = vec![vec![T::zero(); k]; k];
    for p in 0..k {
        for q in 0..k {
            omega_tilde[p][q] = gamma_ratio
                * s.factor_covariance_omega
                    .get(p)
                    .and_then(|r| r.get(q))
                    .copied()
                    .unwrap_or(T::zero());
        }
        let c = s
            .factor_crowding_scores
            .get(p)
            .copied()
            .unwrap_or(T::zero());
        omega_tilde[p][p] += s.crowding_penalty_lambda / gamma_i * c * c;
    }
    let mut mtx = vec![vec![T::zero(); k]; k];
    let mut rhs = vec![T::zero(); k];
    for p in 0..k {
        rhs[p] = (0..k)
            .map(|q| omega_tilde[p][q] * ba[q])
            .fold(T::zero(), |a, b| a + b);
        for q in 0..k {
            mtx[p][q] = if p == q { T::one() } else { T::zero() };
            for r in 0..k {
                mtx[p][q] += omega_tilde[p][r] * g[r][q];
            }
        }
    }
    let beta_tilde = if k == 0 {
        Vec::new()
    } else {
        paleo_solve(mtx, rhs)
    };
    let crowding_penalized_weights: Vec<T> = (0..n)
        .map(|i| {
            let factor = (0..k)
                .map(|p| {
                    s.factor_loadings_b
                        .get(i)
                        .and_then(|r| r.get(p))
                        .copied()
                        .unwrap_or(T::zero())
                        * beta_tilde.get(p).copied().unwrap_or(T::zero())
                })
                .fold(T::zero(), |a, b| a + b);
            (s.raw_alphas[i] - factor) / gamma_i / d[i].max_val(T::from_f64(1e-12))
        })
        .collect();
    let expected_active_return = s
        .raw_alphas
        .iter()
        .zip(&crowding_penalized_weights)
        .map(|(a, w)| *a * *w)
        .fold(T::zero(), |a, b| a + b);
    let idiosyncratic_information_ratio = orthogonal_alphas
        .iter()
        .enumerate()
        .map(|(i, x)| *x * *x / d[i].max_val(T::from_f64(1e-12)))
        .fold(T::zero(), |a, b| a + b)
        .sqrt();
    let mut residual_factor_exposures = vec![T::zero(); k];
    for p in 0..k {
        for i in 0..n {
            residual_factor_exposures[p] += s
                .factor_loadings_b
                .get(i)
                .and_then(|r| r.get(p))
                .copied()
                .unwrap_or(T::zero())
                * crowding_penalized_weights[i];
        }
    }
    let num = crowding_penalized_weights
        .iter()
        .enumerate()
        .map(|(i, w)| *w * d[i] * *w)
        .fold(T::zero(), |a, b| a + b);
    let mut den = num;
    for i in 0..n {
        for j in 0..n {
            let mut bo = T::zero();
            for p in 0..k {
                for q in 0..k {
                    bo += s
                        .factor_loadings_b
                        .get(i)
                        .and_then(|r| r.get(p))
                        .copied()
                        .unwrap_or(T::zero())
                        * s.factor_covariance_omega
                            .get(p)
                            .and_then(|r| r.get(q))
                            .copied()
                            .unwrap_or(T::zero())
                        * s.factor_loadings_b
                            .get(j)
                            .and_then(|r| r.get(q))
                            .copied()
                            .unwrap_or(T::zero());
                }
            }
            den += crowding_penalized_weights[i] * crowding_penalized_weights[j] * bo;
        }
    }
    PaleologoFactorNeutralAlphaReport {
        orthogonal_alphas,
        shrunk_idio_variances: d,
        pure_factor_neutral_weights,
        crowding_penalized_weights,
        residual_factor_exposures,
        idiosyncratic_information_ratio,
        proportion_idiosyncratic_variance_pvi: num / den.max_val(T::from_f64(1e-12)),
        expected_active_return,
    }
}

#[derive(Debug, Error, Clone, Copy, PartialEq, Eq)]
pub enum AllocationError {
    #[error("singular covariance or view matrix")]
    SingularCovariance,
    #[error("invalid matrix dimensions")]
    InvalidDimensions,
    #[error("infeasible box bounds")]
    InfeasibleBounds,
}

#[derive(Clone, Debug)]
pub struct BoxConstrainedMarkowitzSpec<T: Real> {
    pub covariance: Vec<Vec<T>>,
    pub expected_returns: Vec<T>,
    pub risk_aversion: T,
    pub lower_bounds: Vec<T>,
    pub upper_bounds: Vec<T>,
}

pub fn try_box_constrained_markowitz_weights<T: Real>(
    s: &BoxConstrainedMarkowitzSpec<T>,
) -> Result<Vec<T>, AllocationError> {
    let n = s.expected_returns.len();
    if n == 0
        || s.covariance.len() != n
        || s.covariance.iter().any(|r| r.len() != n)
        || s.lower_bounds.len() != n
        || s.upper_bounds.len() != n
        || !s.risk_aversion.val().is_finite()
        || s.covariance
            .iter()
            .flat_map(|r| r.iter())
            .any(|x| !x.val().is_finite())
        || s.expected_returns.iter().any(|x| !x.val().is_finite())
        || s.lower_bounds
            .iter()
            .chain(s.upper_bounds.iter())
            .any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    if s.risk_aversion.val() <= 0.0 {
        return Err(AllocationError::InvalidDimensions);
    }
    if s.lower_bounds
        .iter()
        .zip(&s.upper_bounds)
        .any(|(l, u)| l.val() > u.val() + 1e-12)
        || s.lower_bounds.iter().map(|x| x.val()).sum::<f64>() > 1.0 + 1e-12
        || s.upper_bounds.iter().map(|x| x.val()).sum::<f64>() < 1.0 - 1e-12
    {
        return Err(AllocationError::InfeasibleBounds);
    }
    // Determine the active set on primal f64 data.  No iterative solve is
    // recorded on an AD tape; the reduced KKT system is lifted once below.
    let mut state = vec![0i8; n]; // -1 lower, +1 upper, 0 free
    let cf = s
        .covariance
        .iter()
        .map(|r| r.iter().map(|x| x.val()).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let muf = s
        .expected_returns
        .iter()
        .map(|x| x.val())
        .collect::<Vec<_>>();
    let lf = s.lower_bounds.iter().map(|x| x.val()).collect::<Vec<_>>();
    let uf = s.upper_bounds.iter().map(|x| x.val()).collect::<Vec<_>>();
    let gam = s.risk_aversion.val();
    let solve_primal = |state: &[i8]| -> Result<(Vec<f64>, f64), AllocationError> {
        let free = (0..n).filter(|&i| state[i] == 0).collect::<Vec<_>>();
        if free.is_empty() {
            let w = state
                .iter()
                .enumerate()
                .map(|(i, z)| if *z < 0 { lf[i] } else { uf[i] })
                .collect::<Vec<_>>();
            let g_fn = |i: usize| muf[i] - gam * (0..n).map(|j| cf[i][j] * w[j]).sum::<f64>();
            let lo_lam = (0..n)
                .filter(|&i| state[i] > 0)
                .map(g_fn)
                .fold(f64::NEG_INFINITY, f64::max);
            let hi_lam = (0..n)
                .filter(|&i| state[i] < 0)
                .map(g_fn)
                .fold(f64::INFINITY, f64::min);
            let lambda = match (lo_lam.is_finite(), hi_lam.is_finite()) {
                (true, true) => 0.5 * (lo_lam + hi_lam),
                (true, false) => lo_lam,
                (false, true) => hi_lam,
                (false, false) => 0.0,
            };
            return Ok((w, lambda));
        }
        let m = free.len() + 1;
        let mut a = vec![vec![0.0; m]; m];
        let mut b = vec![0.0; m];
        let fixed_sum = (0..n)
            .filter(|&i| state[i] != 0)
            .map(|i| if state[i] < 0 { lf[i] } else { uf[i] })
            .sum::<f64>();
        for (ii, &i) in free.iter().enumerate() {
            for (jj, &j) in free.iter().enumerate() {
                a[ii][jj] = gam * cf[i][j];
            }
            a[ii][m - 1] = 1.0;
            a[m - 1][ii] = 1.0;
            b[ii] = muf[i]
                - (0..n)
                    .filter(|&j| state[j] != 0)
                    .map(|j| gam * cf[i][j] * if state[j] < 0 { lf[j] } else { uf[j] })
                    .sum::<f64>();
        }
        b[m - 1] = 1.0 - fixed_sum;
        let sol = solve_linear_generic::<f64>(&a, &b)?;
        let mut w = vec![0.0; n];
        for i in 0..n {
            w[i] = if state[i] < 0 {
                lf[i]
            } else if state[i] > 0 {
                uf[i]
            } else {
                0.0
            };
        }
        for (ii, &i) in free.iter().enumerate() {
            w[i] = sol[ii];
        }
        Ok((w, sol[m - 1]))
    };
    for _ in 0..(4 * n + 16) {
        let (w, lambda) = solve_primal(&state)?;
        let fixed_sum = (0..n).filter(|&i| state[i] != 0).map(|i| w[i]).sum::<f64>();
        if state.iter().all(|z| *z != 0) && (fixed_sum - 1.0).abs() > 1e-10 {
            // An all-active partition can only be feasible when the bounds
            // themselves sum to one.  Release one bound and re-solve the
            // equality KKT system instead of returning an infeasible budget.
            let want_upper = fixed_sum > 1.0;
            if let Some(i) =
                (0..n).find(|&i| (want_upper && state[i] > 0) || (!want_upper && state[i] < 0))
            {
                state[i] = 0;
                continue;
            }
        }
        let mut worst: Option<(usize, f64)> = None;
        for i in 0..n {
            if state[i] == 0 {
                let v = if w[i] < lf[i] - 1e-11 {
                    lf[i] - w[i]
                } else if w[i] > uf[i] + 1e-11 {
                    w[i] - uf[i]
                } else {
                    0.0
                };
                if v > 0.0 && worst.as_ref().is_none_or(|(_, x)| v > *x) {
                    worst = Some((i, v));
                }
            }
        }
        if let Some((i, _)) = worst {
            state[i] = if w[i] < lf[i] { -1 } else { 1 };
            continue;
        }
        // Release an incorrectly active bound using the KKT multiplier.
        let mut release: Option<(usize, f64)> = None;
        for i in 0..n {
            if state[i] == 0 {
                continue;
            }
            let grad = gam * (0..n).map(|j| cf[i][j] * w[j]).sum::<f64>() - muf[i] + lambda;
            let violation = if state[i] < 0 {
                (-grad).max(0.0)
            } else {
                grad.max(0.0)
            };
            if violation > 1e-10 && release.is_none_or(|(_, x)| violation > x) {
                release = Some((i, violation));
            }
        }
        if let Some((i, _)) = release {
            state[i] = 0;
            continue;
        }
        break;
    }
    let _ = solve_primal(&state)?;
    let free = (0..n).filter(|&i| state[i] == 0).collect::<Vec<_>>();
    let mut w = vec![T::zero(); n];
    let fixed = (0..n).filter(|&i| state[i] != 0).fold(T::zero(), |a, i| {
        a + if state[i] < 0 {
            s.lower_bounds[i]
        } else {
            s.upper_bounds[i]
        }
    });
    for i in 0..n {
        if state[i] < 0 {
            w[i] = s.lower_bounds[i];
        } else if state[i] > 0 {
            w[i] = s.upper_bounds[i];
        }
    }
    if free.is_empty() {
        return Ok(w);
    }
    if free.len() == 1 {
        w[free[0]] = T::one() - fixed;
        return Ok(w);
    }
    let m = free.len() + 1;
    let mut a = vec![vec![T::zero(); m]; m];
    let mut b = vec![T::zero(); m];
    for (ii, &i) in free.iter().enumerate() {
        for (jj, &j) in free.iter().enumerate() {
            a[ii][jj] = s.risk_aversion * s.covariance[i][j];
        }
        a[ii][m - 1] = T::one();
        a[m - 1][ii] = T::one();
        b[ii] = s.expected_returns[i]
            - (0..n)
                .filter(|&j| state[j] != 0)
                .map(|j| s.risk_aversion * s.covariance[i][j] * w[j])
                .fold(T::zero(), |x, y| x + y);
    }
    b[m - 1] = T::one() - fixed;
    let sol = solve_linear_generic(&a, &b)?;
    for (ii, &i) in free.iter().enumerate() {
        w[i] = sol[ii];
    }
    Ok(w)
}
pub fn box_constrained_markowitz_weights<T: Real>(s: &BoxConstrainedMarkowitzSpec<T>) -> Vec<T> {
    try_box_constrained_markowitz_weights(s).expect("feasible Markowitz bounds")
}

pub fn try_long_only_markowitz_weights<T: Real>(
    covariance: &[Vec<T>],
    expected_returns: &[T],
    risk_aversion: T,
) -> Result<Vec<T>, AllocationError> {
    let n = expected_returns.len();
    try_box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
        covariance: covariance.to_vec(),
        expected_returns: expected_returns.to_vec(),
        risk_aversion,
        lower_bounds: vec![T::zero(); n],
        upper_bounds: vec![T::one(); n],
    })
}

pub fn long_only_markowitz_weights<T: Real>(
    covariance: &[Vec<T>],
    expected_returns: &[T],
    risk_aversion: T,
) -> Vec<T> {
    try_long_only_markowitz_weights(covariance, expected_returns, risk_aversion)
        .expect("feasible long-only Markowitz inputs")
}

#[derive(Clone, Debug)]
pub struct BlackLittermanSpec<T: Real> {
    pub covariance: Vec<Vec<T>>,
    pub market_weights: Vec<T>,
    pub risk_aversion: T,
    pub tau: T,
    pub view_matrix: Vec<Vec<T>>,
    pub view_returns: Vec<T>,
    pub omega: Vec<Vec<T>>,
}
#[derive(Clone, Debug)]
pub struct BlackLittermanPosterior<T: Real> {
    pub expected_returns: Vec<T>,
    pub covariance: Vec<Vec<T>>,
    pub weights: Vec<T>,
}
impl<T: Real> BlackLittermanSpec<T> {
    pub fn new(
        covariance: Vec<Vec<T>>,
        market_weights: Vec<T>,
        risk_aversion: T,
        tau: T,
        view_matrix: Vec<Vec<T>>,
        view_returns: Vec<T>,
        omega: Vec<Vec<T>>,
    ) -> Self {
        Self {
            covariance,
            market_weights,
            risk_aversion,
            tau,
            view_matrix,
            view_returns,
            omega,
        }
    }
}

/// AD-preserving Gaussian elimination with partial pivoting.  Pivot choice is
/// primal-only, but every row operation is performed in `T`, so Dual/Dual2
/// and TapeVar derivatives survive the solve.
pub fn solve_linear_generic<T: Real>(a: &[Vec<T>], b: &[T]) -> Result<Vec<T>, AllocationError> {
    let n = b.len();
    if n == 0
        || a.len() != n
        || a.iter().any(|r| r.len() != n)
        || a.iter().flatten().any(|x| !x.val().is_finite())
        || b.iter().any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let mut m = vec![vec![T::zero(); n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            m[i][j] = a[i][j];
        }
        m[i][n] = b[i];
    }
    for k in 0..n {
        let Some(p) = (k..n).max_by(|&i, &j| m[i][k].val().abs().total_cmp(&m[j][k].val().abs()))
        else {
            return Err(AllocationError::SingularCovariance);
        };
        if !m[p][k].val().is_finite() || m[p][k].val().abs() <= 1e-14 {
            return Err(AllocationError::SingularCovariance);
        }
        m.swap(k, p);
        let pivot = m[k][k];
        for value in m[k].iter_mut().take(n + 1).skip(k) {
            *value /= pivot;
        }
        for i in 0..n {
            if i == k {
                continue;
            }
            let f = m[i][k];
            let pivot_row = m[k].clone();
            for j in k..=n {
                m[i][j] -= f * pivot_row[j];
            }
        }
    }
    Ok(m.into_iter().map(|r| r[n]).collect())
}

#[derive(Clone, Debug)]
pub struct HierarchicalRiskParityResult<T: Real> {
    pub weights: Vec<T>,
    pub seriated_order: Vec<usize>,
}
fn hrp_validate<T: Real>(cov: &[Vec<T>]) -> Result<(), AllocationError> {
    let n = cov.len();
    if n == 0
        || cov.iter().any(|r| r.len() != n)
        || cov.iter().flatten().any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    for i in 0..n {
        if cov[i][i].val() <= 0.0 {
            return Err(AllocationError::SingularCovariance);
        }
        for j in 0..n {
            if (cov[i][j].val() - cov[j][i].val()).abs() > 1e-10 {
                return Err(AllocationError::InvalidDimensions);
            }
            let corr = cov[i][j].val() / (cov[i][i].val() * cov[j][j].val()).sqrt();
            if !corr.is_finite() || corr.abs() > 1.0 + 1e-8 {
                return Err(AllocationError::InvalidDimensions);
            }
        }
    }
    Ok(())
}
pub fn try_hierarchical_risk_parity<T: Real>(
    cov: &[Vec<T>],
) -> Result<HierarchicalRiskParityResult<T>, AllocationError> {
    hrp_validate(cov)?;
    let n = cov.len();
    let mut raw_d = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let corr =
                (cov[i][j].val() / (cov[i][i].val() * cov[j][j].val()).sqrt()).clamp(-1.0, 1.0);
            raw_d[i][j] = (0.5 * (1.0 - corr)).max(0.0).sqrt();
        }
    }
    let mut euclid_d = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let d = (0..n)
                .map(|k| {
                    let x = raw_d[i][k] - raw_d[j][k];
                    x * x
                })
                .sum::<f64>()
                .sqrt();
            euclid_d[i][j] = d;
            euclid_d[j][i] = d;
        }
    }
    let mut clusters: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    let dist = |a: &[usize], b: &[usize]| {
        let mut out = f64::INFINITY;
        for &i in a {
            for &j in b {
                out = out.min(euclid_d[i][j]);
            }
        }
        out
    };
    while clusters.len() > 1 {
        let mut best = (0usize, 1usize, f64::INFINITY);
        for i in 0..clusters.len() {
            for j in i + 1..clusters.len() {
                let d = dist(&clusters[i], &clusters[j]);
                if d < best.2 {
                    best = (i, j, d);
                }
            }
        }
        let right = clusters.remove(best.1);
        let mut left = clusters.remove(best.0);
        left.extend(right);
        clusters.push(left);
    }
    let order = clusters.pop().unwrap_or_default();
    let seriated_order = order.clone();
    let mut weights = vec![T::one(); n];
    let mut active = vec![order];
    while let Some(l) = active.pop() {
        if l.len() < 2 {
            continue;
        }
        let m = l.len() / 2;
        let (a, b) = (l[..m].to_vec(), l[m..].to_vec());
        let variance = |c: &[usize]| {
            let mut inv = T::zero();
            for &i in c {
                inv += T::one() / cov[i][i];
            }
            let mut v = T::zero();
            for &i in c {
                for &j in c {
                    v += (T::one() / cov[i][i] / inv) * cov[i][j] * (T::one() / cov[j][j] / inv);
                }
            }
            v
        };
        let va = variance(&a);
        let vb = variance(&b);
        let den = (va + vb).max_val(T::from_f64(1e-14));
        for &i in &a {
            weights[i] *= vb / den;
        }
        for &i in &b {
            weights[i] *= va / den;
        }
        active.push(a);
        active.push(b);
    }
    let sum = weights.iter().copied().fold(T::zero(), |x, y| x + y);
    for w in &mut weights {
        *w /= sum.max_val(T::from_f64(1e-14));
    }
    Ok(HierarchicalRiskParityResult {
        weights,
        seriated_order,
    })
}
pub fn hierarchical_risk_parity_weights<T: Real>(cov: &[Vec<T>]) -> Vec<T> {
    try_hierarchical_risk_parity(cov)
        .expect("valid HRP covariance")
        .weights
}

#[derive(Clone, Debug)]
pub struct EntropyPoolingSpec<T: Real = f64> {
    pub scenario_returns: Vec<Vec<T>>,
    pub prior_probs: Vec<T>,
    pub view_matrix: Vec<Vec<T>>,
    pub view_targets: Vec<T>,
    pub risk_aversion: T,
}

#[derive(Clone, Debug)]
pub struct EntropyPoolingResult<T: Real> {
    pub posterior_probs: Vec<T>,
    pub posterior_mean: Vec<T>,
    pub posterior_covariance: Vec<Vec<T>>,
    pub weights: Vec<T>,
    pub kl_divergence: T,
    pub effective_num_scenarios: T,
}

pub fn try_entropy_pooling_posterior<T: Real>(
    s: &EntropyPoolingSpec<T>,
) -> Result<EntropyPoolingResult<T>, AllocationError> {
    let n = s.prior_probs.len();
    let k = s.view_matrix.len();
    if n == 0
        || s.scenario_returns.len() != n
        || s.scenario_returns[0].is_empty()
        || s.scenario_returns.iter().any(|r| {
            r.len() != s.scenario_returns[0].len() || r.iter().any(|x| !x.val().is_finite())
        })
        || s.prior_probs
            .iter()
            .any(|x| !x.val().is_finite() || x.val() <= 0.0)
        || s.view_matrix
            .iter()
            .any(|r| r.len() != n || r.iter().any(|x| !x.val().is_finite()))
        || s.view_targets.len() != k
        || s.view_targets.iter().any(|x| !x.val().is_finite())
        || !s.risk_aversion.val().is_finite()
        || s.risk_aversion.val() <= 0.0
    {
        return Err(AllocationError::InvalidDimensions);
    }

    // 1. Primal f64 Newton solve for Lagrange multipliers \lambda^*
    let mut lam_f64 = vec![0.0_f64; k];
    if k > 0 {
        for _ in 0..2 {
            let mut logits = vec![0.0_f64; n];
            for i in 0..n {
                for j in 0..k {
                    logits[i] += lam_f64[j] * s.view_matrix[j][i].val();
                }
            }
            let mx = logits.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let mut p_f64 = vec![0.0_f64; n];
            let mut z = 0.0_f64;
            for i in 0..n {
                p_f64[i] = s.prior_probs[i].val() * (logits[i] - mx).exp();
                z += p_f64[i];
            }
            for x in &mut p_f64 {
                *x /= z.max(1e-16);
            }
            let mut g = vec![0.0_f64; k];
            for j in 0..k {
                for i in 0..n {
                    g[j] += p_f64[i] * s.view_matrix[j][i].val();
                }
                g[j] -= s.view_targets[j].val();
            }
            if g.iter().all(|x| x.abs() < 1e-13) {
                break;
            }
            let mut h = vec![vec![0.0_f64; k]; k];
            let mut mu = vec![0.0_f64; k];
            for j in 0..k {
                for i in 0..n {
                    mu[j] += p_f64[i] * s.view_matrix[j][i].val();
                }
            }
            for a in 0..k {
                for b in 0..k {
                    for i in 0..n {
                        h[a][b] += p_f64[i]
                            * (s.view_matrix[a][i].val() - mu[a])
                            * (s.view_matrix[b][i].val() - mu[b]);
                    }
                    if a == b {
                        h[a][b] += 1e-10;
                    }
                }
            }
            let rhs: Vec<f64> = g.iter().map(|&x| -x).collect();
            let step = solve_linear_generic(&h, &rhs)?;
            for j in 0..k {
                lam_f64[j] += step[j];
            }
        }
    }

    // 2. Lift \lambda^* into T: Real via 2 generic Newton IFT polish steps
    let mut lam: Vec<T> = lam_f64.into_iter().map(T::from_f64).collect();
    let mut p = vec![T::zero(); n];
    for _ in 0..12 {
        let mut logits = vec![T::zero(); n];
        for i in 0..n {
            for j in 0..k {
                logits[i] += lam[j] * s.view_matrix[j][i];
            }
        }
        let mx = logits
            .iter()
            .map(|x| x.val())
            .fold(f64::NEG_INFINITY, f64::max);
        let mut z = T::zero();
        for i in 0..n {
            p[i] = s.prior_probs[i] * (logits[i] - mx).exp();
            z += p[i];
        }
        for x in &mut p {
            *x /= z.max_val(T::from_f64(1e-16));
        }
        if k > 0 {
            let mut g = vec![T::zero(); k];
            let mut mu = vec![T::zero(); k];
            for j in 0..k {
                for i in 0..n {
                    mu[j] += p[i] * s.view_matrix[j][i];
                }
                g[j] = mu[j] - s.view_targets[j];
            }
            let mut h = vec![vec![T::zero(); k]; k];
            for a in 0..k {
                for b in 0..k {
                    for i in 0..n {
                        h[a][b] +=
                            p[i] * (s.view_matrix[a][i] - mu[a]) * (s.view_matrix[b][i] - mu[b]);
                    }
                    if a == b {
                        h[a][b] += T::from_f64(1e-10);
                    }
                }
            }
            let rhs: Vec<T> = g.iter().map(|&x| -x).collect();
            let step = solve_linear_generic(&h, &rhs)?;
            for j in 0..k {
                lam[j] += step[j];
            }
        }
    }

    let m = s.scenario_returns[0].len();
    let mut mean = vec![T::zero(); m];
    for i in 0..n {
        for j in 0..m {
            mean[j] += p[i] * s.scenario_returns[i][j];
        }
    }
    let mut cov = vec![vec![T::zero(); m]; m];
    let mut reg_cov = vec![vec![T::zero(); m]; m];
    for i in 0..n {
        for a in 0..m {
            for b in 0..m {
                let term = p[i]
                    * (s.scenario_returns[i][a] - mean[a])
                    * (s.scenario_returns[i][b] - mean[b]);
                cov[a][b] += term;
                reg_cov[a][b] += term * s.risk_aversion;
            }
        }
    }
    for (a, row) in reg_cov.iter_mut().enumerate().take(m) {
        row[a] += T::from_f64(1e-6);
    }
    let weights = solve_linear_generic(&reg_cov, &mean)?;

    let kl = p
        .iter()
        .enumerate()
        .map(|(i, &x)| {
            x * (x.max_val(T::from_f64(1e-16)) / s.prior_probs[i].max_val(T::from_f64(1e-16))).ln()
        })
        .fold(T::zero(), |a, b| a + b);
    let en = p
        .iter()
        .map(|&x| {
            let q = x.max_val(T::from_f64(1e-16));
            -q * q.ln()
        })
        .fold(T::zero(), |a, b| a + b)
        .exp();

    Ok(EntropyPoolingResult {
        posterior_probs: p,
        posterior_mean: mean,
        posterior_covariance: cov,
        weights,
        kl_divergence: kl,
        effective_num_scenarios: en,
    })
}

pub fn entropy_pooling_posterior<T: Real>(s: &EntropyPoolingSpec<T>) -> EntropyPoolingResult<T> {
    try_entropy_pooling_posterior(s).expect("invalid entropy pooling inputs")
}

#[derive(Clone, Debug)]
pub struct SchurHrpSpec<T: Real = f64> {
    pub covariance: Vec<Vec<T>>,
    pub gamma: T,
}

#[derive(Clone, Debug)]
pub struct SchurHrpReport<T: Real> {
    pub weights: Vec<T>,
    pub portfolio_variance: T,
}

fn invert_spd_generic<T: Real>(mat: &[Vec<T>]) -> Result<Vec<Vec<T>>, AllocationError> {
    let n = mat.len();
    let mut cols = vec![vec![T::zero(); n]; n];
    for j in 0..n {
        let mut e = vec![T::zero(); n];
        e[j] = T::one();
        let col = solve_linear_generic(mat, &e)?;
        for i in 0..n {
            cols[i][j] = col[i];
        }
    }
    Ok(cols)
}

fn schur_bisect_weights<T: Real>(
    cov: &[Vec<T>],
    rhs: &[T],
    gamma: T,
) -> Result<Vec<T>, AllocationError> {
    let n = cov.len();
    if n == 0 || rhs.len() != n {
        return Err(AllocationError::InvalidDimensions);
    }
    if n == 1 {
        return Ok(vec![T::one()]);
    }
    let n1 = n / 2;
    let n2 = n - n1;
    let mut a = vec![vec![T::zero(); n1]; n1];
    let mut b = vec![vec![T::zero(); n2]; n1];
    let mut c = vec![vec![T::zero(); n1]; n2];
    let mut d = vec![vec![T::zero(); n2]; n2];
    for i in 0..n1 {
        for j in 0..n1 {
            a[i][j] = cov[i][j];
        }
        for j in 0..n2 {
            b[i][j] = cov[i][n1 + j];
            c[j][i] = cov[n1 + j][i];
        }
    }
    for i in 0..n2 {
        for j in 0..n2 {
            d[i][j] = cov[n1 + i][n1 + j];
        }
    }

    let a_inv = invert_spd_generic(&a)?;
    let d_inv = invert_spd_generic(&d)?;

    let mut bd_inv = vec![vec![T::zero(); n2]; n1];
    for i in 0..n1 {
        for j in 0..n2 {
            for k in 0..n2 {
                bd_inv[i][j] += b[i][k] * d_inv[k][j];
            }
        }
    }
    let mut a_c = a;
    let mut rhs_a = vec![T::zero(); n1];
    for i in 0..n1 {
        let mut bd_r = T::zero();
        for k in 0..n2 {
            bd_r += bd_inv[i][k] * rhs[n1 + k];
        }
        rhs_a[i] = (rhs[i] - gamma * bd_r).max_val(T::from_f64(1e-4));
        for j in 0..n1 {
            let mut bdc_ij = T::zero();
            for k in 0..n2 {
                bdc_ij += bd_inv[i][k] * c[k][j];
            }
            a_c[i][j] -= gamma * bdc_ij;
        }
    }

    // Compute C A^{-1} B and C A^{-1} 1_{n1}
    let mut ca_inv = vec![vec![T::zero(); n1]; n2];
    for i in 0..n2 {
        for j in 0..n1 {
            for k in 0..n1 {
                ca_inv[i][j] += c[i][k] * a_inv[k][j];
            }
        }
    }
    let mut d_c = d;
    let mut rhs_d = vec![T::zero(); n2];
    for i in 0..n2 {
        let mut ca_r = T::zero();
        for k in 0..n1 {
            ca_r += ca_inv[i][k] * rhs[k];
        }
        rhs_d[i] = (rhs[n1 + i] - gamma * ca_r).max_val(T::from_f64(1e-4));
        for j in 0..n2 {
            let mut cab_ij = T::zero();
            for k in 0..n1 {
                cab_ij += ca_inv[i][k] * b[k][j];
            }
            d_c[i][j] -= gamma * cab_ij;
        }
    }

    let w_a = schur_bisect_weights(&a_c, &rhs_a, gamma)?;
    let w_d = schur_bisect_weights(&d_c, &rhs_d, gamma)?;
    let mut quad_a = T::zero();
    let mut dot_ra = T::zero();
    for i in 0..n1 {
        dot_ra += w_a[i] * rhs_a[i];
        for j in 0..n1 {
            quad_a += w_a[i] * a_c[i][j] * w_a[j];
        }
    }
    let mut quad_d = T::zero();
    let mut dot_rd = T::zero();
    for i in 0..n2 {
        dot_rd += w_d[i] * rhs_d[i];
        for j in 0..n2 {
            quad_d += w_d[i] * d_c[i][j] * w_d[j];
        }
    }
    let nu_a = quad_a / dot_ra.max_val(T::from_f64(1e-16));
    let nu_d = quad_d / dot_rd.max_val(T::from_f64(1e-16));
    let alpha_a = nu_d / (nu_a + nu_d).max_val(T::from_f64(1e-16));
    let mut out = Vec::with_capacity(n);
    for x in w_a {
        out.push(alpha_a * x);
    }
    for x in w_d {
        out.push((T::one() - alpha_a) * x);
    }
    Ok(out)
}

pub fn try_schur_complementary_hrp<T: Real>(
    spec: &SchurHrpSpec<T>,
) -> Result<SchurHrpReport<T>, AllocationError> {
    let n = spec.covariance.len();
    if n == 0 || spec.covariance.iter().any(|r| r.len() != n) {
        return Err(AllocationError::InvalidDimensions);
    }
    let g = spec.gamma.max_val(T::zero()).min_val(T::one());
    let ones = vec![T::one(); n];
    let weights = schur_bisect_weights(&spec.covariance, &ones, g)?;
    let mut port_var = T::zero();
    for i in 0..n {
        for j in 0..n {
            port_var += weights[i] * spec.covariance[i][j] * weights[j];
        }
    }
    Ok(SchurHrpReport {
        weights,
        portfolio_variance: port_var,
    })
}

pub fn schur_complementary_hrp<T: Real>(spec: &SchurHrpSpec<T>) -> SchurHrpReport<T> {
    try_schur_complementary_hrp(spec).expect("valid Schur-HRP covariance")
}

#[derive(Clone, Debug)]
pub struct OasNonlinearShrinkageReport<T: Real> {
    pub oas_covariance: Vec<Vec<T>>,
    pub oas_shrinkage_intensity: T,
    pub nonlinear_shrunk_covariance: Vec<Vec<T>>,
}

pub fn try_oas_and_nonlinear_shrinkage_covariance<T: Real>(
    returns: &[Vec<T>],
) -> Result<OasNonlinearShrinkageReport<T>, AllocationError> {
    let n_obs = returns.len();
    if n_obs < 2 {
        return Err(AllocationError::InvalidDimensions);
    }
    let p = returns[0].len();
    if p == 0 || returns.iter().any(|r| r.len() != p) {
        return Err(AllocationError::InvalidDimensions);
    }
    let mut mean = vec![T::zero(); p];
    for row in returns {
        for j in 0..p {
            mean[j] += row[j] / n_obs as f64;
        }
    }
    let mut s = vec![vec![T::zero(); p]; p];
    let denom = (n_obs - 1) as f64;
    for row in returns {
        for i in 0..p {
            let di = row[i] - mean[i];
            for j in 0..p {
                s[i][j] += di * (row[j] - mean[j]) / denom;
            }
        }
    }

    let mut tr_s = T::zero();
    let mut tr_s2 = T::zero();
    for i in 0..p {
        tr_s += s[i][i];
        for j in 0..p {
            tr_s2 += s[i][j] * s[j][i];
        }
    }
    let p_f = p as f64;
    let n_f = n_obs as f64;
    let num = tr_s2 * (1.0 - 2.0 / p_f) + tr_s * tr_s;
    let den = (tr_s2 - tr_s * tr_s / p_f).max_val(T::from_f64(1e-14)) * (n_f + 1.0 - 2.0 / p_f);
    let rho_oas = (num / den).max_val(T::zero()).min_val(T::one());
    let mu_target = tr_s / p_f;

    let mut oas_cov = vec![vec![T::zero(); p]; p];
    for i in 0..p {
        for j in 0..p {
            oas_cov[i][j] = s[i][j] * (T::one() - rho_oas);
            if i == j {
                oas_cov[i][j] += rho_oas * mu_target;
            }
        }
    }

    // Ledoit-Wolf (2012/2020) nonlinear spectral resolvent shrinkage:
    // S_NL = (1 - c_nl) * S + c_nl * mu_target * (S / mu_target + I)^{-1} * (1 + mu_target)
    let c = (p_f / n_f).min(0.95);
    let mut reg = vec![vec![T::zero(); p]; p];
    for i in 0..p {
        for j in 0..p {
            reg[i][j] = s[i][j] / mu_target.max_val(T::from_f64(1e-12));
            if i == j {
                reg[i][j] += T::one();
            }
        }
    }
    let reg_inv = invert_spd_generic(&reg)?;
    let mut nl_raw = vec![vec![T::zero(); p]; p];
    let mut tr_nl = T::zero();
    for i in 0..p {
        for j in 0..p {
            let diag_term = if i == j {
                mu_target * (2.0 * c)
            } else {
                T::zero()
            };
            nl_raw[i][j] =
                s[i][j] * (1.0 - c) + diag_term - reg_inv[i][j] * mu_target * (c * (2.0 - c));
        }
        tr_nl += nl_raw[i][i];
    }
    let scale_tr = tr_s / tr_nl.max_val(T::from_f64(1e-14));
    let mut nl_cov = vec![vec![T::zero(); p]; p];
    for i in 0..p {
        for j in 0..p {
            nl_cov[i][j] = (nl_raw[i][j] + nl_raw[j][i]) * 0.5 * scale_tr;
        }
    }

    Ok(OasNonlinearShrinkageReport {
        oas_covariance: oas_cov,
        oas_shrinkage_intensity: rho_oas,
        nonlinear_shrunk_covariance: nl_cov,
    })
}

pub fn oas_and_nonlinear_shrinkage_covariance<T: Real>(
    returns: &[Vec<T>],
) -> OasNonlinearShrinkageReport<T> {
    try_oas_and_nonlinear_shrinkage_covariance(returns).expect("valid return matrix for OAS")
}

pub fn implied_equilibrium_returns<T: Real>(
    cov: &[Vec<T>],
    weights: &[T],
    risk_aversion: T,
) -> Vec<T> {
    (0..cov.len())
        .map(|i| {
            risk_aversion
                * cov[i]
                    .iter()
                    .zip(weights)
                    .map(|(a, b)| *a * *b)
                    .fold(T::zero(), |a, b| a + b)
        })
        .collect()
}
pub fn try_black_litterman_posterior<T: Real>(
    s: &BlackLittermanSpec<T>,
) -> Result<BlackLittermanPosterior<T>, AllocationError> {
    let n = s.covariance.len();
    let k = s.view_returns.len();
    if n == 0
        || s.market_weights.len() != n
        || s.covariance.iter().any(|r| r.len() != n)
        || s.view_matrix.len() != k
        || s.view_matrix.iter().any(|r| r.len() != n)
        || s.omega.len() != k
        || s.omega.iter().any(|r| r.len() != k)
        || !s.risk_aversion.val().is_finite()
        || s.risk_aversion.val() <= 0.0
        || !s.tau.val().is_finite()
        || s.tau.val() < 0.0
        || s.market_weights.iter().any(|x| !x.val().is_finite())
        || s.covariance
            .iter()
            .flat_map(|r| r.iter())
            .any(|x| !x.val().is_finite())
        || s.view_matrix
            .iter()
            .flat_map(|r| r.iter())
            .chain(s.view_returns.iter())
            .chain(s.omega.iter().flat_map(|r| r.iter()))
            .any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let pi = implied_equilibrium_returns(&s.covariance, &s.market_weights, s.risk_aversion);
    if k == 0 {
        return Ok(BlackLittermanPosterior {
            expected_returns: pi.clone(),
            covariance: s.covariance.clone(),
            weights: s.market_weights.clone(),
        });
    }
    let tau_cov = s
        .covariance
        .iter()
        .map(|row| row.iter().map(|x| *x * s.tau).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut middle = vec![vec![T::zero(); k]; k];
    for i in 0..k {
        for j in 0..k {
            middle[i][j] = s.omega[i][j]
                + (0..n)
                    .flat_map(|a| (0..n).map(move |b| (a, b)))
                    .map(|(a, b)| s.view_matrix[i][a] * tau_cov[a][b] * s.view_matrix[j][b])
                    .fold(T::zero(), |x, y| x + y);
        }
    }
    let mut rhs_y = vec![vec![T::zero(); n]; k];
    for i in 0..k {
        for j in 0..n {
            rhs_y[i][j] = (0..n)
                .map(|a| s.view_matrix[i][a] * tau_cov[a][j])
                .fold(T::zero(), |x, y| x + y);
        }
    }
    let mut ymat = vec![vec![T::zero(); n]; k];
    for j in 0..n {
        let rhs = (0..k).map(|i| rhs_y[i][j]).collect::<Vec<_>>();
        let sol = solve_linear_generic(&middle, &rhs)?;
        for i in 0..k {
            ymat[i][j] = sol[i];
        }
    }
    let rhs_v = (0..k)
        .map(|i| {
            s.view_returns[i]
                - (0..n)
                    .map(|j| s.view_matrix[i][j] * pi[j])
                    .fold(T::zero(), |x, y| x + y)
        })
        .collect::<Vec<_>>();
    let v = solve_linear_generic(&middle, &rhs_v)?;
    let mu = (0..n)
        .map(|i| {
            pi[i]
                + (0..k)
                    .map(|a| {
                        (0..n)
                            .map(|j| tau_cov[i][j] * s.view_matrix[a][j])
                            .fold(T::zero(), |x, y| x + y)
                            * v[a]
                    })
                    .fold(T::zero(), |x, y| x + y)
        })
        .collect::<Vec<_>>();
    let mut post_cov = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            post_cov[i][j] = s.covariance[i][j] + tau_cov[i][j]
                - (0..k)
                    .map(|a| {
                        (0..n)
                            .map(|b| tau_cov[i][b] * s.view_matrix[a][b])
                            .fold(T::zero(), |x, y| x + y)
                            * ymat[a][j]
                    })
                    .fold(T::zero(), |x, y| x + y);
        }
    }
    let weighted = s
        .covariance
        .iter()
        .map(|row| row.iter().map(|x| *x * s.risk_aversion).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let raw = solve_linear_generic(&weighted, &mu)?;
    let z = raw.iter().fold(T::zero(), |a, x| a + *x);
    if !z.val().is_finite() || z.val().abs() <= 1e-14 {
        return Err(AllocationError::SingularCovariance);
    }
    let w = raw.iter().map(|x| *x / z).collect::<Vec<_>>();
    Ok(BlackLittermanPosterior {
        expected_returns: mu,
        covariance: post_cov,
        weights: w,
    })
}
pub fn black_litterman_posterior<T: Real>(s: &BlackLittermanSpec<T>) -> BlackLittermanPosterior<T> {
    try_black_litterman_posterior(s).expect("valid Black-Litterman covariance/view system")
}
pub fn try_markowitz_minimum_variance_weights<T: Real>(
    cov: &[Vec<T>],
) -> Result<Vec<T>, AllocationError> {
    let n = cov.len();
    if n == 0 || cov.iter().any(|r| r.len() != n) {
        return Err(AllocationError::InvalidDimensions);
    }
    let ones = vec![T::one(); n];
    let x = solve_linear_generic(cov, &ones)?;
    let z = x.iter().fold(T::zero(), |a, b| a + *b);
    if !z.val().is_finite() || z.val().abs() <= 1e-14 {
        return Err(AllocationError::SingularCovariance);
    }
    Ok(x.into_iter().map(|x| x / z).collect())
}
pub fn markowitz_minimum_variance_weights<T: Real>(cov: &[Vec<T>]) -> Vec<T> {
    try_markowitz_minimum_variance_weights(cov).expect("valid covariance matrix")
}
pub fn try_markowitz_tangency_weights<T: Real>(
    cov: &[Vec<T>],
    excess: &[T],
) -> Result<Vec<T>, AllocationError> {
    let n = cov.len();
    if n == 0 || excess.len() != n || cov.iter().any(|r| r.len() != n) {
        return Err(AllocationError::InvalidDimensions);
    }
    let x = solve_linear_generic(cov, excess)?;
    let z = x.iter().fold(T::zero(), |a, b| a + *b);
    if !z.val().is_finite() || z.val().abs() <= 1e-14 {
        return Err(AllocationError::SingularCovariance);
    }
    Ok(x.into_iter().map(|x| x / z).collect())
}
pub fn markowitz_tangency_weights<T: Real>(cov: &[Vec<T>], excess: &[T]) -> Vec<T> {
    try_markowitz_tangency_weights(cov, excess).expect("valid covariance matrix")
}

fn ccd(c: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut y = b.iter().map(|x| x.sqrt()).collect::<Vec<_>>();
    for _ in 0..5000 {
        let mut max: f64 = 0.0;
        for i in 0..n {
            let off = (0..n)
                .filter(|&j| j != i)
                .map(|j| c[i][j] * y[j])
                .sum::<f64>();
            let disc = off * off + 4. * c[i][i].max(1e-14) * b[i];
            let ny = (-off + disc.sqrt()) / (2. * c[i][i].max(1e-14));
            max = max.max((ny - y[i]).abs());
            y[i] = ny;
        }
        if max < 1e-14 {
            break;
        }
    }
    for _ in 0..5 {
        let f = (0..n)
            .map(|i| (0..n).map(|j| c[i][j] * y[j]).sum::<f64>() - b[i] / y[i])
            .collect::<Vec<_>>();
        let h = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            c[i][j] + b[i] / (y[i] * y[i])
                        } else {
                            c[i][j]
                        }
                    })
                    .collect()
            })
            .collect::<Vec<Vec<_>>>();
        let d = solve_linear_generic::<f64>(&h, &f).expect("risk-budgeting Hessian is SPD");
        for i in 0..n {
            y[i] -= d[i];
        }
    }
    y
}
pub fn try_risk_budgeting_weights<T: Real>(
    cov: &[Vec<T>],
    budgets: &[T],
) -> Result<Vec<T>, AllocationError> {
    let n = cov.len();
    if n == 0
        || budgets.len() != n
        || cov.iter().any(|r| r.len() != n)
        || budgets
            .iter()
            .any(|x| !x.val().is_finite() || x.val() <= 0.0)
        || cov
            .iter()
            .flat_map(|r| r.iter())
            .any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    // Risk budgeting is defined for a symmetric positive-definite covariance
    // matrix.  Without this guard the coordinate descent below can generate
    // NaNs (negative diagonal or an asymmetric input) and those NaNs used to
    // leak out as apparently valid weights.
    for i in 0..n {
        if cov[i][i].val() <= 0.0 {
            return Err(AllocationError::SingularCovariance);
        }
        for j in (i + 1)..n {
            if (cov[i][j].val() - cov[j][i].val()).abs() > 1e-10 {
                return Err(AllocationError::InvalidDimensions);
            }
        }
    }
    let sum = budgets.iter().map(|x| x.val()).sum::<f64>().max(1e-14);
    let b = budgets.iter().map(|x| x.val() / sum).collect::<Vec<_>>();
    let cf = cov
        .iter()
        .map(|r| r.iter().map(|x| x.val()).collect())
        .collect::<Vec<Vec<_>>>();
    let y = ccd(&cf, &b);
    if y.iter().any(|x| !x.is_finite() || *x <= 0.0) {
        return Err(AllocationError::SingularCovariance);
    }
    let h = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        cf[i][j] + b[i] / (y[i] * y[i])
                    } else {
                        cf[i][j]
                    }
                })
                .collect()
        })
        .collect::<Vec<Vec<_>>>();
    let inv = (0..n)
        .map(|i| {
            solve_linear_generic::<f64>(
                &h,
                &(0..n)
                    .map(|j| if i == j { 1. } else { 0. })
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let budget_sum_t = budgets.iter().fold(T::zero(), |a, x| a + *x);
    let mut yt = vec![T::zero(); n];
    for i in 0..n {
        let d = (0..n)
            .map(|j| {
                T::from_f64(-inv[j][i])
                    * (((0..n)
                        .map(|k| cov[j][k] * y[k])
                        .fold(T::zero(), |a, x| a + x)
                        - budgets[j] / budget_sum_t / y[j])
                        - T::from_f64((0..n).map(|k| cf[j][k] * y[k]).sum::<f64>() - b[j] / y[j]))
            })
            .fold(T::zero(), |a, x| a + x);
        yt[i] = T::from_f64(y[i]) + d;
    }
    let z = yt.iter().fold(T::zero(), |a, x| a + *x);
    if !z.val().is_finite() || z.val().abs() <= 1e-14 {
        return Err(AllocationError::SingularCovariance);
    }
    Ok(yt.into_iter().map(|x| x / z).collect())
}
pub fn risk_budgeting_weights<T: Real>(cov: &[Vec<T>], budgets: &[T]) -> Vec<T> {
    try_risk_budgeting_weights(cov, budgets).expect("valid risk-budgeting covariance and budgets")
}
pub fn equal_risk_contribution_weights<T: Real>(cov: &[Vec<T>]) -> Vec<T> {
    let n = cov.len();
    if n == 0 {
        return Vec::new();
    }
    risk_budgeting_weights(cov, &vec![T::one() / n as f64; n])
}
pub fn euler_risk_contributions<T: Real>(cov: &[Vec<T>], w: &[T]) -> (T, Vec<T>) {
    let mut q = T::zero();
    for i in 0..w.len() {
        for j in 0..w.len() {
            q += w[i] * cov[i][j] * w[j];
        }
    }
    let sig = q.sqrt();
    let rc = (0..w.len())
        .map(|i| {
            w[i] * (0..w.len())
                .map(|j| cov[i][j] * w[j])
                .fold(T::zero(), |a, x| a + x)
                / sig
        })
        .collect();
    (sig, rc)
}

#[derive(Clone, Debug)]
pub struct MeanCvarSpec<T: Real> {
    pub scenario_returns: Vec<Vec<T>>,
    pub expected_returns: Vec<T>,
    pub risk_aversion: T,
    pub confidence: f64,
    pub smoothing_eps: f64,
    pub lower_bounds: Vec<T>,
    pub upper_bounds: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct MeanCvarResult<T: Real> {
    pub weights: Vec<T>,
    pub var_threshold: T,
    pub cvar: T,
    pub expected_return: T,
}
fn project_capped_simplex<T: Real>(
    v: &[T],
    lower: &[T],
    upper: &[T],
) -> Result<Vec<T>, AllocationError> {
    let n = v.len();
    if n == 0
        || lower.len() != n
        || upper.len() != n
        || lower
            .iter()
            .zip(upper)
            .any(|(l, u)| l.val() > u.val() + 1e-12)
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let lo_sum = lower.iter().copied().fold(T::zero(), |a, z| a + z);
    let hi_sum = upper.iter().copied().fold(T::zero(), |a, z| a + z);
    if lo_sum.val() > 1.0 + 1e-12 || hi_sum.val() < 1.0 - 1e-12 {
        return Err(AllocationError::InfeasibleBounds);
    }
    let mut lo = v
        .iter()
        .zip(upper)
        .map(|(x, u)| x.val() - u.val())
        .fold(f64::INFINITY, f64::min);
    let mut hi = v
        .iter()
        .zip(lower)
        .map(|(x, l)| x.val() - l.val())
        .fold(f64::NEG_INFINITY, f64::max);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let sum: f64 = v
            .iter()
            .zip(lower)
            .zip(upper)
            .map(|((x, l), u)| (x.val() - mid).max(l.val()).min(u.val()))
            .sum();
        if sum > 1.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let mid = 0.5 * (lo + hi);
    let mut free_sum = T::zero();
    let mut clamped_sum = T::zero();
    let mut free_cnt = 0usize;
    for i in 0..v.len() {
        let cand = v[i].val() - mid;
        if cand <= lower[i].val() {
            clamped_sum += lower[i];
        } else if cand >= upper[i].val() {
            clamped_sum += upper[i];
        } else {
            free_sum += v[i];
            free_cnt += 1;
        }
    }
    let shift = if free_cnt > 0 {
        (free_sum + clamped_sum - T::one()) / (free_cnt as f64)
    } else {
        T::from_f64(mid)
    };
    Ok((0..v.len())
        .map(|i| {
            let cand = v[i].val() - mid;
            if cand <= lower[i].val() {
                lower[i]
            } else if cand >= upper[i].val() {
                upper[i]
            } else {
                (v[i] - shift).max_val(lower[i]).min_val(upper[i])
            }
        })
        .collect())
}

pub fn try_mean_cvar_weights<T: Real>(
    s: &MeanCvarSpec<T>,
) -> Result<MeanCvarResult<T>, AllocationError> {
    let n = s.expected_returns.len();
    let m_scen = s.scenario_returns.len();
    if n == 0
        || m_scen == 0
        || s.scenario_returns.iter().any(|r| r.len() != n)
        || s.lower_bounds.len() != n
        || s.upper_bounds.len() != n
        || !(0.5..1.0).contains(&s.confidence)
        || !s.risk_aversion.val().is_finite()
        || s.risk_aversion.val() < 0.0
        || !s.smoothing_eps.is_finite()
        || s.smoothing_eps < 0.0
        || s.scenario_returns
            .iter()
            .flat_map(|r| r.iter())
            .chain(s.expected_returns.iter())
            .chain(s.lower_bounds.iter())
            .chain(s.upper_bounds.iter())
            .any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let eps = T::from_f64(s.smoothing_eps.max(1e-5));
    let tail_mass = T::from_f64(((1.0 - s.confidence) * m_scen as f64).max(1e-12));
    let idx = ((s.confidence * m_scen as f64).floor() as usize).min(m_scen - 1);
    let init_guess: Vec<T> = (0..n).map(|_| T::from_f64(1.0 / n as f64)).collect();
    let mut w = project_capped_simplex(&init_guess, &s.lower_bounds, &s.upper_bounds)?;
    let solve_zeta = |weights: &[T]| -> (T, Vec<T>) {
        let losses: Vec<T> = s
            .scenario_returns
            .iter()
            .map(|r| {
                -r.iter()
                    .zip(weights)
                    .map(|(x, y)| *x * *y)
                    .fold(T::zero(), |a, z| a + z)
            })
            .collect();
        let mut sorted = losses.clone();
        sorted.sort_by(|a, b| a.val().total_cmp(&b.val()));
        let mut zeta = sorted[idx];
        for _ in 0..4 {
            let mut f = -tail_mass;
            let mut df = T::zero();
            for &l in &losses {
                let u = l - zeta;
                let rad = (u * u + eps * eps).sqrt();
                f += (T::one() + u / rad) * 0.5;
                df += eps * eps / (rad * rad * rad) * 0.5;
            }
            zeta += f / df.max_val(T::from_f64(1e-12));
        }
        (zeta, losses)
    };
    let eval_obj = |weights: &[T]| -> T {
        let (z, ls) = solve_zeta(weights);
        let mut h_sum = T::zero();
        for &l in &ls {
            let u = l - z;
            h_sum += (u + (u * u + eps * eps).sqrt()) * 0.5;
        }
        let cv = z + h_sum / tail_mass;
        let mu = s
            .expected_returns
            .iter()
            .zip(weights)
            .map(|(a, b)| *a * *b)
            .fold(T::zero(), |a, x| a + x);
        -mu + s.risk_aversion * cv
    };
    for _ in 0..32 {
        let (zeta, losses) = solve_zeta(&w);
        let mut grad = vec![T::zero(); n];
        let mut g_norm2 = 0.0;
        for i in 0..n {
            let mut tail_g = T::zero();
            for m in 0..m_scen {
                let u = losses[m] - zeta;
                let dphi = (T::one() + u / (u * u + eps * eps).sqrt()) * 0.5;
                tail_g -= (dphi / tail_mass) * s.scenario_returns[m][i];
            }
            grad[i] = -s.expected_returns[i] + s.risk_aversion * tail_g;
            g_norm2 += grad[i].val() * grad[i].val();
        }
        let g_norm = g_norm2.sqrt().max(1e-12);
        let base_obj = eval_obj(&w).val();
        let mut step = 0.5 / g_norm;
        for _ in 0..6 {
            let trial_in: Vec<T> = (0..n).map(|i| w[i] - grad[i] * step).collect();
            let cand = project_capped_simplex(&trial_in, &s.lower_bounds, &s.upper_bounds)?;
            if eval_obj(&cand).val() <= base_obj + 1e-12 {
                w = cand;
                break;
            }
            step *= 0.5;
        }
    }
    let (zeta, losses) = solve_zeta(&w);
    let mut smooth_hinge_sum = T::zero();
    for &l in &losses {
        let u = l - zeta;
        smooth_hinge_sum += (u + (u * u + eps * eps).sqrt()) * 0.5;
    }
    let cvar = zeta + smooth_hinge_sum / tail_mass;
    let er = s
        .expected_returns
        .iter()
        .zip(&w)
        .map(|(a, b)| *a * *b)
        .fold(T::zero(), |a, z| a + z);
    Ok(MeanCvarResult {
        weights: w,
        var_threshold: zeta,
        cvar,
        expected_return: er,
    })
}

pub fn mean_cvar_weights<T: Real>(s: &MeanCvarSpec<T>) -> MeanCvarResult<T> {
    try_mean_cvar_weights(s).expect("valid mean-CVaR inputs")
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum ShrinkageTarget {
    ScaledIdentity,
    ConstantVarianceIdentity,
    DiagonalVariance,
    ConstantCorrelation,
}
#[derive(Clone, Debug)]
pub struct LedoitWolfResult<T: Real> {
    pub shrunk_covariance: Vec<Vec<T>>,
    pub sample_covariance: Vec<Vec<T>>,
    pub target_matrix: Vec<Vec<T>>,
    pub shrinkage_intensity: T,
}
pub fn ledoit_wolf_shrinkage<T: Real>(
    r: &[Vec<T>],
    target: ShrinkageTarget,
) -> Result<LedoitWolfResult<T>, AllocationError> {
    let n = r.len();
    if n < 2
        || r.iter().any(|x| {
            x.is_empty() || x.len() != r[0].len() || x.iter().any(|v| !v.val().is_finite())
        })
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let p = r[0].len();
    let mean = (0..p)
        .map(|j| r.iter().map(|x| x[j]).fold(T::zero(), |a, z| a + z) / n as f64)
        .collect::<Vec<_>>();
    let s = (0..p)
        .map(|i| {
            (0..p)
                .map(|j| {
                    r.iter()
                        .map(|x| (x[i] - mean[i]) * (x[j] - mean[j]))
                        .fold(T::zero(), |a, z| a + z)
                        / n as f64
                })
                .collect()
        })
        .collect::<Vec<Vec<_>>>();
    let mut f = vec![vec![T::zero(); p]; p];
    match target {
        ShrinkageTarget::ScaledIdentity | ShrinkageTarget::ConstantVarianceIdentity => {
            let m = (0..p).map(|i| s[i][i]).fold(T::zero(), |a, z| a + z) / p as f64;
            for (i, row) in f.iter_mut().enumerate().take(p) {
                row[i] = m;
            }
        }
        ShrinkageTarget::DiagonalVariance => {
            for (i, row) in f.iter_mut().enumerate().take(p) {
                row[i] = s[i][i];
            }
        }
        ShrinkageTarget::ConstantCorrelation => {
            let mut c = T::zero();
            let mut cnt: f64 = 0.;
            for i in 0..p {
                for j in 0..i {
                    c += s[i][j] / (s[i][i] * s[j][j]).sqrt().max_val(T::from_f64(1e-12));
                    cnt += 1.;
                }
            }
            let c = c / cnt.max(1.);
            for i in 0..p {
                for j in 0..p {
                    f[i][j] = if i == j {
                        s[i][i]
                    } else {
                        c * (s[i][i] * s[j][j]).sqrt()
                    };
                }
            }
        }
    }
    // Ledoit--Wolf intensity: estimate the variance of the sample
    // covariance elementwise (pi-hat) and divide by the squared distance to
    // the target (gamma).  The former fixed `0.1/gamma` heuristic was not a
    // shrinkage estimator and could materially over/under-shrink depending
    // only on units.
    let mut pi = T::zero();
    for row in r {
        for i in 0..p {
            for j in 0..p {
                if i == j
                    && matches!(
                        target,
                        ShrinkageTarget::DiagonalVariance | ShrinkageTarget::ConstantCorrelation
                    )
                {
                    continue;
                }
                let outer = (row[i] - mean[i]) * (row[j] - mean[j]);
                let err = outer - s[i][j];
                pi += err * err;
            }
        }
    }
    pi = pi / n as f64;
    let d = (0..p)
        .map(|i| {
            (0..p)
                .map(|j| (s[i][j] - f[i][j]) * (s[i][j] - f[i][j]))
                .fold(T::zero(), |a, z| a + z)
        })
        .fold(T::zero(), |a, z| a + z);
    let delta = ((pi / n as f64) / (d + T::from_f64(1e-12)))
        .max_val(T::zero())
        .min_val(T::one());
    let shr = (0..p)
        .map(|i| {
            (0..p)
                .map(|j| delta * f[i][j] + (T::one() - delta) * s[i][j])
                .collect()
        })
        .collect();
    Ok(LedoitWolfResult {
        shrunk_covariance: shr,
        sample_covariance: s,
        target_matrix: f,
        shrinkage_intensity: delta,
    })
}
#[derive(Clone, Debug)]
pub struct RobustMarkowitzSpec<T: Real> {
    pub covariance: Vec<Vec<T>>,
    pub return_uncertainty_cov: Vec<Vec<T>>,
    pub expected_returns: Vec<T>,
    pub risk_aversion: T,
    pub uncertainty_radius: T,
    pub lower_bounds: Vec<T>,
    pub upper_bounds: Vec<T>,
}
pub fn try_robust_markowitz_weights<T: Real>(
    s: &RobustMarkowitzSpec<T>,
) -> Result<Vec<T>, AllocationError> {
    let n = s.expected_returns.len();
    // Validate every matrix/vector before the fixed-point loop.  The old
    // implementation only checked `n`, so malformed uncertainty matrices
    // could panic while indexing (and negative radii produced a nonsensical
    // adversarial covariance adjustment).
    if n == 0
        || s.covariance.len() != n
        || s.covariance.iter().any(|r| r.len() != n)
        || s.return_uncertainty_cov.len() != n
        || s.return_uncertainty_cov.iter().any(|r| r.len() != n)
        || s.lower_bounds.len() != n
        || s.upper_bounds.len() != n
        || !s.risk_aversion.val().is_finite()
        || s.risk_aversion.val() <= 0.0
        || !s.uncertainty_radius.val().is_finite()
        || s.uncertainty_radius.val() < 0.0
        || s.expected_returns
            .iter()
            .chain(s.lower_bounds.iter())
            .chain(s.upper_bounds.iter())
            .any(|x| !x.val().is_finite())
        || s.covariance
            .iter()
            .chain(s.return_uncertainty_cov.iter())
            .flat_map(|r| r.iter())
            .any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let mut cov = s.covariance.clone();
    let mut w = try_box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
        covariance: cov.clone(),
        expected_returns: s.expected_returns.clone(),
        risk_aversion: s.risk_aversion,
        lower_bounds: s.lower_bounds.clone(),
        upper_bounds: s.upper_bounds.clone(),
    })?;
    for iter in 0..32 {
        let mut ow = T::zero();
        for i in 0..n {
            for j in 0..n {
                ow += w[i] * s.return_uncertainty_cov[i][j] * w[j];
            }
        }
        let scale = s.uncertainty_radius / (ow + T::from_f64(1e-16)).sqrt();
        for i in 0..n {
            for j in 0..n {
                cov[i][j] = s.covariance[i][j]
                    + scale * s.return_uncertainty_cov[i][j]
                        / s.risk_aversion.max_val(T::from_f64(1e-12));
            }
        }
        let nw = try_box_constrained_markowitz_weights(&BoxConstrainedMarkowitzSpec {
            covariance: cov.clone(),
            expected_returns: s.expected_returns.clone(),
            risk_aversion: s.risk_aversion,
            lower_bounds: s.lower_bounds.clone(),
            upper_bounds: s.upper_bounds.clone(),
        })?;
        if iter >= 16
            && nw
                .iter()
                .zip(&w)
                .all(|(a, b)| (a.val() - b.val()).abs() < 1e-14)
        {
            w = nw;
            break;
        }
        w = nw;
    }
    Ok(w)
}
pub fn robust_markowitz_weights<T: Real>(s: &RobustMarkowitzSpec<T>) -> Vec<T> {
    try_robust_markowitz_weights(s).expect("valid robust Markowitz inputs")
}

pub fn gerber_covariance_matrix<T: Real>(
    returns: &[Vec<T>],
    threshold_c: f64,
    smoothing_eps: f64,
) -> Result<Vec<Vec<T>>, AllocationError> {
    let m = returns.len();
    if m < 2
        || returns[0].is_empty()
        || returns
            .iter()
            .any(|r| r.len() != returns[0].len() || r.iter().any(|x| !x.val().is_finite()))
        || !(0.0..2.0).contains(&threshold_c)
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let n = returns[0].len();
    let mut mean = vec![T::zero(); n];
    for r in returns {
        for j in 0..n {
            mean[j] += r[j] / m as f64;
        }
    }
    let mut sd = vec![T::zero(); n];
    for j in 0..n {
        let mut q = T::zero();
        for r in returns {
            let d = r[j] - mean[j];
            q += d * d;
        }
        sd[j] = (q / (m as f64 - 1.0)).sqrt().max_val(T::from_f64(1e-12));
    }
    let eps = smoothing_eps.max(1e-4);
    let step = |u: T| (T::one() + u / (u * u + eps * eps).sqrt()) * 0.5;
    let mut out = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        out[i][i] = sd[i] * sd[i];
        for j in i + 1..n {
            let (mut con, mut dis) = (T::zero(), T::zero());
            for r in returns {
                let zi = (r[i] - mean[i]) / sd[i];
                let zj = (r[j] - mean[j]) / sd[j];
                let ui = step(zi - threshold_c);
                let di = step(-zi - threshold_c);
                let uj = step(zj - threshold_c);
                let dj = step(-zj - threshold_c);
                con += ui * uj + di * dj;
                dis += ui * dj + di * uj;
            }
            let g = (con - dis) / (con + dis + 1e-12);
            out[i][j] = sd[i] * sd[j] * g;
            out[j][i] = out[i][j];
        }
    }
    Ok(out)
}

#[allow(clippy::needless_range_loop)]
pub fn marchenko_pastur_denoise_covariance<T: Real>(
    cov: &[Vec<T>],
    num_observations: usize,
) -> Result<Vec<Vec<T>>, AllocationError> {
    let n = cov.len();
    if n == 0
        || num_observations < 2
        || cov.iter().any(|r| r.len() != n)
        || cov.iter().flatten().any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    for i in 0..n {
        if cov[i][i].val() <= 0.0 {
            return Err(AllocationError::SingularCovariance);
        }
        for j in i + 1..n {
            if (cov[i][j].val() - cov[j][i].val()).abs() > 1e-10 {
                return Err(AllocationError::SingularCovariance);
            }
        }
    }
    let sd: Vec<T> = (0..n).map(|i| cov[i][i].sqrt()).collect();
    let mut c = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            c[i][j] = cov[i][j] / (sd[i] * sd[j]);
        }
    }
    let mut v = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        v[i][i] = T::one();
    }
    for _ in 0..(32 * n * n).max(64) {
        let mut p = (0, 1);
        let mut mx = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                if c[i][j].val().abs() > mx {
                    mx = c[i][j].val().abs();
                    p = (i, j);
                }
            }
        }
        if mx < 1e-14 {
            break;
        }
        let (i, j) = p;
        let tau = (c[j][j] - c[i][i]) / (c[i][j] * 2.0);
        let sgn = if tau.val() >= 0.0 { 1.0 } else { -1.0 };
        let t = T::from_f64(sgn) / (tau.abs() + (T::one() + tau * tau).sqrt());
        let co = T::one() / (T::one() + t * t).sqrt();
        let si = t * co;
        for k in 0..n {
            let a = c[k][i];
            let b = c[k][j];
            c[k][i] = co * a - si * b;
            c[k][j] = si * a + co * b;
        }
        for k in 0..n {
            let a = c[i][k];
            let b = c[j][k];
            c[i][k] = co * a - si * b;
            c[j][k] = si * a + co * b;
        }
        for k in 0..n {
            let a = v[k][i];
            let b = v[k][j];
            v[k][i] = co * a - si * b;
            v[k][j] = si * a + co * b;
        }
    }
    let mut eig: Vec<T> = (0..n).map(|i| c[i][i]).collect();
    let lp = (1.0 + (n as f64 / num_observations as f64).sqrt()).powi(2);
    let mut ns = T::zero();
    let mut cnt = 0;
    for x in &eig {
        if x.val() <= lp {
            ns += *x;
            cnt += 1;
        }
    }
    if cnt > 0 {
        let m = ns / cnt as f64;
        for x in &mut eig {
            if x.val() <= lp {
                *x = m;
            }
        }
    }
    let mut out = vec![vec![T::zero(); n]; n];
    let mut raw = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                raw[i][j] += eig[k] * (v[i][k] * v[j][k]);
            }
        }
    }
    let ds: Vec<T> = (0..n)
        .map(|i| raw[i][i].max_val(T::from_f64(1e-14)).sqrt())
        .collect();
    for i in 0..n {
        for j in 0..n {
            out[i][j] = sd[i] * sd[j] * raw[i][j] / (ds[i] * ds[j]);
        }
        out[i][i] = cov[i][i];
    }
    Ok(out)
}
#[derive(Clone, Debug)]
pub struct NestedClusteredOptimizationSpec<T: Real = f64> {
    pub covariance: Vec<Vec<T>>,
    pub expected_returns: Option<Vec<T>>,
    pub risk_aversion: T,
    pub num_clusters: usize,
}
#[allow(clippy::needless_range_loop)]
pub fn try_nested_clustered_optimization_weights<T: Real>(
    s: &NestedClusteredOptimizationSpec<T>,
) -> Result<Vec<T>, AllocationError> {
    let n = s.covariance.len();
    if n == 0
        || s.covariance.iter().any(|r| r.len() != n)
        || s.num_clusters == 0
        || s.num_clusters > n
        || s.covariance.iter().flatten().any(|x| !x.val().is_finite())
    {
        return Err(AllocationError::InvalidDimensions);
    }
    if let Some(mu) = &s.expected_returns
        && (mu.len() != n || mu.iter().any(|x| !x.val().is_finite()))
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let mut sd = vec![0.0; n];
    for i in 0..n {
        if s.covariance[i][i].val() <= 0.0 {
            return Err(AllocationError::SingularCovariance);
        }
        sd[i] = s.covariance[i][i].val().sqrt();
        for j in i + 1..n {
            if (s.covariance[i][j].val() - s.covariance[j][i].val()).abs() > 1e-10 {
                return Err(AllocationError::SingularCovariance);
            }
        }
    }
    let k = s.num_clusters;
    let mut clusters: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
    while clusters.len() > k {
        let (mut ba, mut bb, mut bd) = (0, 1, f64::INFINITY);
        for a in 0..clusters.len() {
            for b in a + 1..clusters.len() {
                let mut d = f64::INFINITY;
                for &i in &clusters[a] {
                    for &j in &clusters[b] {
                        d = d.min(
                            (0.5 * (1.0
                                - (s.covariance[i][j].val() / (sd[i] * sd[j])).clamp(-1.0, 1.0)))
                            .max(0.0)
                            .sqrt(),
                        );
                    }
                }
                if d < bd {
                    (ba, bb, bd) = (a, b, d);
                }
            }
        }
        let right = clusters.remove(bb);
        clusters[ba].extend(right);
    }
    let gamma = s.risk_aversion.max_val(T::from_f64(1e-12));
    let mut local = vec![vec![T::zero(); n]; k];
    for (ci, cl) in clusters.iter().enumerate() {
        let m = cl.len();
        let mut sub = vec![vec![T::zero(); m]; m];
        let mut rhs = vec![T::one(); m];
        for (a, &i) in cl.iter().enumerate() {
            if let Some(mu) = &s.expected_returns {
                rhs[a] = T::one() + mu[i] / gamma;
            }
            for (b, &j) in cl.iter().enumerate() {
                sub[a][b] = s.covariance[i][j];
            }
        }
        let x = solve_linear_generic(&sub, &rhs)?;
        let z = x.iter().fold(T::zero(), |a, b| a + *b);
        if z.val().abs() < 1e-14 {
            return Err(AllocationError::SingularCovariance);
        }
        for (a, &i) in cl.iter().enumerate() {
            local[ci][i] = x[a] / z;
        }
    }
    let mut red = vec![vec![T::zero(); k]; k];
    let mut rhs = vec![T::one(); k];
    for a in 0..k {
        if let Some(mu) = &s.expected_returns {
            let mut cm = T::zero();
            for &i in &clusters[a] {
                cm += local[a][i] * mu[i];
            }
            rhs[a] = T::one() + cm / gamma;
        }
        for b in 0..k {
            for &i in &clusters[a] {
                for &j in &clusters[b] {
                    red[a][b] += local[a][i] * s.covariance[i][j] * local[b][j];
                }
            }
        }
    }
    let om = solve_linear_generic(&red, &rhs)?;
    let z = om.iter().fold(T::zero(), |a, b| a + *b);
    if z.val().abs() < 1e-14 {
        return Err(AllocationError::SingularCovariance);
    }
    let mut out = vec![T::zero(); n];
    for (ci, cl) in clusters.iter().enumerate() {
        for &i in cl {
            out[i] = om[ci] / z * local[ci][i];
        }
    }
    Ok(out)
}
pub fn nested_clustered_optimization_weights<T: Real>(
    s: &NestedClusteredOptimizationSpec<T>,
) -> Vec<T> {
    try_nested_clustered_optimization_weights(s).expect("valid NCO inputs")
}

#[derive(Clone, Debug)]
pub struct IdzorekBlackLittermanSpec<T: Real = f64> {
    pub covariance: Vec<Vec<T>>,
    pub market_weights: Vec<T>,
    pub risk_aversion: T,
    pub tau: T,
    pub view_matrix: Vec<Vec<T>>,
    pub view_returns: Vec<T>,
    pub view_confidences: Vec<T>,
}

pub fn try_idzorek_black_litterman_posterior<T: Real>(
    s: &IdzorekBlackLittermanSpec<T>,
) -> Result<BlackLittermanPosterior<T>, AllocationError> {
    let n = s.covariance.len();
    let k = s.view_returns.len();
    if n == 0
        || s.market_weights.len() != n
        || s.market_weights.iter().any(|x| !x.val().is_finite())
        || s.covariance
            .iter()
            .any(|r| r.len() != n || r.iter().any(|x| !x.val().is_finite()))
        || s.view_matrix.len() != k
        || s.view_matrix
            .iter()
            .any(|r| r.len() != n || r.iter().any(|x| !x.val().is_finite()))
        || s.view_returns.iter().any(|x| !x.val().is_finite())
        || s.view_confidences.len() != k
        || s.view_confidences.iter().any(|x| !x.val().is_finite())
        || !s.risk_aversion.val().is_finite()
        || s.risk_aversion.val() <= 0.0
        || !s.tau.val().is_finite()
        || s.tau.val() <= 0.0
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let pi = implied_equilibrium_returns(&s.covariance, &s.market_weights, s.risk_aversion);
    if k == 0 {
        return Ok(BlackLittermanPosterior {
            expected_returns: pi,
            covariance: s.covariance.clone(),
            weights: s.market_weights.clone(),
        });
    }

    let mut p_sigma = vec![vec![T::zero(); n]; k];
    for i in 0..k {
        for j in 0..n {
            for a in 0..n {
                p_sigma[i][j] += s.view_matrix[i][a] * s.covariance[a][j];
            }
        }
    }
    let mut p_sigma_pt = vec![vec![T::zero(); k]; k];
    for i in 0..k {
        for j in 0..k {
            for a in 0..n {
                p_sigma_pt[i][j] += p_sigma[i][a] * s.view_matrix[j][a];
            }
        }
    }

    let mut m_c = vec![vec![T::zero(); k]; k];
    let mut rhs_v = vec![T::zero(); k];
    for i in 0..k {
        let c_i = s.view_confidences[i].max_val(T::zero()).min_val(T::one());
        let one_minus_c = T::one() - c_i;
        for j in 0..k {
            m_c[i][j] = c_i * p_sigma_pt[i][j];
        }
        m_c[i][i] += one_minus_c * p_sigma_pt[i][i].max_val(T::from_f64(1e-14));
        let mut p_pi = T::zero();
        for a in 0..n {
            p_pi += s.view_matrix[i][a] * pi[a];
        }
        rhs_v[i] = c_i * (s.view_returns[i] - p_pi);
    }

    let v = solve_linear_generic(&m_c, &rhs_v)?;
    let mut mu = pi;
    for a in 0..n {
        for i in 0..k {
            mu[a] += p_sigma[i][a] * v[i];
        }
    }

    let mut y_mat = vec![vec![T::zero(); n]; k];
    for j in 0..n {
        let col_rhs: Vec<T> = (0..k)
            .map(|i| {
                let c_i = s.view_confidences[i].max_val(T::zero()).min_val(T::one());
                c_i * p_sigma[i][j]
            })
            .collect();
        let sol = solve_linear_generic(&m_c, &col_rhs)?;
        for i in 0..k {
            y_mat[i][j] = sol[i];
        }
    }

    let mut post_cov = vec![vec![T::zero(); n]; n];
    for a in 0..n {
        for b in 0..n {
            let mut red = T::zero();
            for i in 0..k {
                red += p_sigma[i][a] * y_mat[i][b];
            }
            post_cov[a][b] = s.covariance[a][b] * (T::one() + s.tau) - s.tau * red;
        }
    }

    let raw_w = solve_linear_generic(&s.covariance, &mu)?;
    let weights: Vec<T> = raw_w.into_iter().map(|w| w / s.risk_aversion).collect();
    Ok(BlackLittermanPosterior {
        expected_returns: mu,
        covariance: post_cov,
        weights,
    })
}

pub fn idzorek_black_litterman_posterior<T: Real>(
    s: &IdzorekBlackLittermanSpec<T>,
) -> BlackLittermanPosterior<T> {
    try_idzorek_black_litterman_posterior(s).expect("valid Idzorek Black-Litterman inputs")
}

#[derive(Clone, Debug)]
pub struct KellyAllocationSpec<T: Real = f64> {
    pub expected_returns: Vec<T>,
    pub covariance: Vec<Vec<T>>,
    pub risk_free_rate: T,
    pub max_kelly_fraction: T,
    pub max_drawdown_limit: T,
    pub drawdown_tail_prob: T,
}

#[derive(Clone, Debug)]
pub struct KellyAllocationResult<T: Real> {
    pub weights: Vec<T>,
    pub cash_weight: T,
    pub effective_kelly_fraction: T,
    pub expected_log_growth_rate: T,
    pub portfolio_volatility: T,
}

pub fn try_kelly_fractional_allocation<T: Real>(
    s: &KellyAllocationSpec<T>,
) -> Result<KellyAllocationResult<T>, AllocationError> {
    let n = s.expected_returns.len();
    if n == 0
        || s.expected_returns.iter().any(|x| !x.val().is_finite())
        || s.covariance.len() != n
        || s.covariance
            .iter()
            .any(|r| r.len() != n || r.iter().any(|x| !x.val().is_finite()))
        || !s.risk_free_rate.val().is_finite()
        || !s.max_kelly_fraction.val().is_finite()
        || !s.max_drawdown_limit.val().is_finite()
        || !s.drawdown_tail_prob.val().is_finite()
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let excess: Vec<T> = s
        .expected_returns
        .iter()
        .map(|&m| m - s.risk_free_rate)
        .collect();
    let w_full = solve_linear_generic(&s.covariance, &excess)?;
    let mut sr2 = T::zero();
    for i in 0..n {
        sr2 += excess[i] * w_full[i];
    }
    let sr2_safe = sr2.max_val(T::from_f64(1e-12));
    let d_max = s
        .max_drawdown_limit
        .max_val(T::from_f64(0.01))
        .min_val(T::from_f64(0.95));
    let beta = s
        .drawdown_tail_prob
        .max_val(T::from_f64(1e-4))
        .min_val(T::from_f64(0.99));
    let neg_ln_one_minus_d = -(T::one() - d_max).ln();
    let neg_ln_beta = -beta.ln();
    let f_dd =
        (neg_ln_one_minus_d * 2.0) / (neg_ln_beta + neg_ln_one_minus_d).max_val(T::from_f64(1e-12));
    let f_eff = s.max_kelly_fraction.min_val(f_dd).max_val(T::zero());
    let weights: Vec<T> = w_full.into_iter().map(|w| w * f_eff).collect();
    let sum_w = weights.iter().copied().fold(T::zero(), |a, b| a + b);
    let cash_weight = T::one() - sum_w;
    let growth = s.risk_free_rate + f_eff * sr2_safe - f_eff * f_eff * sr2_safe * 0.5;
    let vol = f_eff * sr2_safe.sqrt();
    Ok(KellyAllocationResult {
        weights,
        cash_weight,
        effective_kelly_fraction: f_eff,
        expected_log_growth_rate: growth,
        portfolio_volatility: vol,
    })
}

pub fn kelly_fractional_allocation<T: Real>(
    s: &KellyAllocationSpec<T>,
) -> KellyAllocationResult<T> {
    try_kelly_fractional_allocation(s).expect("valid Kelly allocation inputs")
}

#[derive(Clone, Debug)]
pub struct FactorEigenRiskParitySpec<T: Real = f64> {
    pub covariance: Vec<Vec<T>>,
    pub factor_risk_budgets: Vec<T>,
}

#[derive(Clone, Debug)]
pub struct FactorEigenRiskParityReport<T: Real = f64> {
    pub weights: Vec<T>,
    pub eigenvalues: Vec<T>,
    pub factor_exposures: Vec<T>,
    pub factor_variance_shares: Vec<T>,
    pub effective_number_of_pca_bets: T,
    pub equal_weight_pca_enb: T,
    pub portfolio_volatility: T,
}

fn jacobi_eigen_spd<T: Real>(a_in: &[Vec<T>]) -> (Vec<T>, Vec<Vec<T>>) {
    let n = a_in.len();
    let mut a = a_in.to_vec();
    let mut v = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        v[i][i] = T::one();
    }
    for _sweep in 0..24 {
        let mut off_norm = 0.0_f64;
        for p in 0..n {
            for q in (p + 1)..n {
                off_norm += a[p][q].val() * a[p][q].val();
            }
        }
        if off_norm.sqrt() < 1e-14 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p][q];
                if apq.val().abs() < 1e-15 {
                    continue;
                }
                let diff = a[q][q] - a[p][p];
                let tau = diff / (apq * 2.0);
                let sign_tau = if tau.val() >= 0.0 {
                    T::one()
                } else {
                    -T::one()
                };
                let t = sign_tau / (tau.abs() + (T::one() + tau * tau).sqrt());
                let c = T::one() / (T::one() + t * t).sqrt();
                let s = t * c;

                let app = a[p][p] - t * apq;
                let aqq = a[q][q] + t * apq;
                a[p][p] = app;
                a[q][q] = aqq;
                a[p][q] = T::zero();
                a[q][p] = T::zero();

                for k in 0..n {
                    if k != p && k != q {
                        let akp = a[k][p];
                        let akq = a[k][q];
                        a[k][p] = c * akp - s * akq;
                        a[p][k] = a[k][p];
                        a[k][q] = s * akp + c * akq;
                        a[q][k] = a[k][q];
                    }
                    let vkp = v[k][p];
                    let vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[j][j].val().total_cmp(&a[i][i].val()));
    let mut evals = Vec::with_capacity(n);
    let mut evecs = vec![vec![T::zero(); n]; n];
    for (col_out, &col_in) in order.iter().enumerate() {
        evals.push(a[col_in][col_in].max_val(T::from_f64(1e-12)));
        let mut col_sum = T::zero();
        for row in 0..n {
            col_sum += v[row][col_in];
        }
        let orient = if col_sum.val() >= 0.0 {
            T::one()
        } else {
            -T::one()
        };
        for row in 0..n {
            evecs[row][col_out] = v[row][col_in] * orient;
        }
    }
    (evals, evecs)
}

/// Roncalli (2013) / Lohre-Opfer-Schiereck (2012) Principal Component Risk Parity
/// (Eigen-Portfolio Factor Risk Budgeting) with Meucci (2009) PCA Diversification Entropy.
pub fn try_solve_principal_component_risk_parity<T: Real>(
    spec: &FactorEigenRiskParitySpec<T>,
) -> Result<FactorEigenRiskParityReport<T>, AllocationError> {
    let n = spec.covariance.len();
    if n == 0
        || spec.factor_risk_budgets.len() != n
        || spec.covariance.iter().any(|r| r.len() != n)
        || spec
            .factor_risk_budgets
            .iter()
            .any(|b| !b.val().is_finite() || b.val() <= 0.0)
    {
        return Err(AllocationError::InvalidDimensions);
    }
    let b_sum = spec
        .factor_risk_budgets
        .iter()
        .fold(T::zero(), |acc, &b| acc + b)
        .max_val(T::from_f64(1e-12));
    let budgets: Vec<T> = spec
        .factor_risk_budgets
        .iter()
        .map(|&b| b / b_sum)
        .collect();

    let (evals, evecs) = jacobi_eigen_spd(&spec.covariance);

    // Closed-form eigen-factor risk budgeting exposures w_tilde_k = sqrt(b_k / lambda_k)
    let w_tilde_raw: Vec<T> = (0..n)
        .map(|k| (budgets[k] / evals[k]).max_val(T::from_f64(1e-14)).sqrt())
        .collect();
    let mut w_raw = vec![T::zero(); n];
    for i in 0..n {
        let mut s = T::zero();
        for k in 0..n {
            s += evecs[i][k] * w_tilde_raw[k];
        }
        w_raw[i] = s;
    }
    let raw_sum = w_raw.iter().fold(T::zero(), |acc, &x| acc + x);
    let weights: Vec<T> = if raw_sum.val().abs() > 1e-8 {
        w_raw.into_iter().map(|x| x / raw_sum).collect()
    } else {
        vec![T::from_f64(1.0 / n as f64); n]
    };

    let factor_exposures: Vec<T> = (0..n)
        .map(|k| (0..n).fold(T::zero(), |acc, i| acc + evecs[i][k] * weights[i]))
        .collect();
    let factor_vars: Vec<T> = (0..n)
        .map(|k| (evals[k] * factor_exposures[k] * factor_exposures[k]).max_val(T::from_f64(1e-16)))
        .collect();
    let tot_var = factor_vars
        .iter()
        .fold(T::zero(), |acc, &v| acc + v)
        .max_val(T::from_f64(1e-14));
    let factor_variance_shares: Vec<T> = factor_vars.iter().map(|&v| v / tot_var).collect();
    let entropy = factor_variance_shares.iter().fold(T::zero(), |acc, &p| {
        acc - p * p.max_val(T::from_f64(1e-14)).ln()
    });
    let enb = entropy.exp();

    // Equal-weight portfolio PCA ENB for comparison
    let ew = T::from_f64(1.0 / n as f64);
    let ew_vars: Vec<T> = (0..n)
        .map(|k| {
            let exp_k = (0..n).fold(T::zero(), |acc, i| acc + evecs[i][k] * ew);
            (evals[k] * exp_k * exp_k).max_val(T::from_f64(1e-16))
        })
        .collect();
    let ew_tot = ew_vars
        .iter()
        .fold(T::zero(), |acc, &v| acc + v)
        .max_val(T::from_f64(1e-14));
    let ew_enb = ew_vars
        .iter()
        .map(|&v| v / ew_tot)
        .fold(T::zero(), |acc, p| {
            acc - p * p.max_val(T::from_f64(1e-14)).ln()
        })
        .exp();

    Ok(FactorEigenRiskParityReport {
        weights,
        eigenvalues: evals,
        factor_exposures,
        factor_variance_shares,
        effective_number_of_pca_bets: enb,
        equal_weight_pca_enb: ew_enb,
        portfolio_volatility: tot_var.sqrt(),
    })
}

pub fn solve_principal_component_risk_parity<T: Real>(
    spec: &FactorEigenRiskParitySpec<T>,
) -> FactorEigenRiskParityReport<T> {
    try_solve_principal_component_risk_parity(spec)
        .expect("valid covariance matrix and factor risk budgets for PCRP")
}

#[derive(Clone, Debug)]
pub struct PottersBouchaudRieCovarianceSpec<T: Real = f64> {
    pub sample_covariance: Vec<Vec<T>>,
    pub num_observations_t: usize,
}
#[derive(Clone, Debug)]
pub struct PottersBouchaudRieCovarianceReport<T: Real = f64> {
    pub aspect_ratio_q: T,
    pub raw_eigenvalues: Vec<T>,
    pub rie_cleaned_eigenvalues: Vec<T>,
    pub raw_condition_number: T,
    pub rie_condition_number: T,
    pub cleaned_covariance: Vec<Vec<T>>,
    pub rie_min_variance_weights: Vec<T>,
    pub raw_min_variance_weights: Vec<T>,
    pub rie_portfolio_variance: T,
}
fn mvw<T: Real>(cov: &[Vec<T>]) -> Vec<T> {
    let inv = invm(cov);
    let n = cov.len();
    let mut z = vec![T::zero(); n];
    for i in 0..n {
        for j in 0..n {
            z[i] += inv[i][j];
        }
    }
    let d = z.iter().copied().fold(T::zero(), |a, x| a + x);
    z.into_iter()
        .map(|x| x / d.max_val(T::from_f64(1e-12)))
        .collect()
}
pub fn solve_potters_bouchaud_rie_nonlinear_shrinkage<T: Real>(
    s: &PottersBouchaudRieCovarianceSpec<T>,
) -> PottersBouchaudRieCovarianceReport<T> {
    let n = s.sample_covariance.len();
    if n == 0 {
        return PottersBouchaudRieCovarianceReport {
            aspect_ratio_q: T::zero(),
            raw_eigenvalues: vec![],
            rie_cleaned_eigenvalues: vec![],
            raw_condition_number: T::zero(),
            rie_condition_number: T::zero(),
            cleaned_covariance: vec![],
            rie_min_variance_weights: vec![],
            raw_min_variance_weights: vec![],
            rie_portfolio_variance: T::zero(),
        };
    }
    let (desc_vals, desc_vecs) = jacobi_eigen_spd(&s.sample_covariance);
    let lambda = desc_vals.into_iter().rev().collect::<Vec<_>>();
    let mut u = vec![vec![T::zero(); n]; n];
    for col in 0..n {
        for row in 0..n {
            u[row][col] = desc_vecs[row][n - 1 - col];
        }
    }
    let q = T::from_f64(n as f64 / s.num_observations_t.max(1) as f64);
    let sqrt_n = T::from_f64(n as f64).sqrt();
    let mut tilde = vec![T::zero(); n];
    for k in 0..n {
        let eta = lambda[k] / sqrt_n;
        let mut gr = T::zero();
        let mut gi = T::zero();
        for j in 0..n {
            let diff = lambda[k] - lambda[j];
            let den = diff * diff + eta * eta;
            gi += eta / den.max_val(T::from_f64(1e-18));
            if j != k {
                gr += diff / den.max_val(T::from_f64(1e-18));
            }
        }
        gr /= T::from_f64(n as f64);
        gi /= T::from_f64(n as f64);
        let re = T::one() - q + q * (lambda[k] * gr + eta * gi);
        let im = q * (lambda[k] * gi - eta * gr);
        tilde[k] = lambda[k] / (re * re + im * im).max_val(T::from_f64(1e-18));
    }
    let trace = lambda.iter().copied().fold(T::zero(), |a, x| a + x);
    let tsum = tilde
        .iter()
        .copied()
        .fold(T::zero(), |a, x| a + x)
        .max_val(T::from_f64(1e-18));
    let clean = tilde
        .iter()
        .map(|x| (*x * trace / tsum).max_val(T::from_f64(1e-12)))
        .collect::<Vec<_>>();
    let mut cov = vec![vec![T::zero(); n]; n];
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                cov[i][j] += clean[k] * u[i][k] * u[j][k];
            }
        }
    }
    let raw_min = lambda
        .iter()
        .copied()
        .fold(T::from_f64(1e12), |a, x| a.min_val(x));
    let raw_max = lambda.iter().copied().fold(T::zero(), |a, x| a.max_val(x));
    let cmin = clean
        .iter()
        .copied()
        .fold(T::from_f64(1e12), |a, x| a.min_val(x));
    let cmax = clean.iter().copied().fold(T::zero(), |a, x| a.max_val(x));
    let w = mvw(&cov);
    let mut pv = T::zero();
    for i in 0..n {
        for j in 0..n {
            pv += w[i] * cov[i][j] * w[j];
        }
    }
    PottersBouchaudRieCovarianceReport {
        aspect_ratio_q: q,
        raw_eigenvalues: lambda,
        rie_cleaned_eigenvalues: clean,
        raw_condition_number: raw_max / raw_min.max_val(T::from_f64(1e-12)),
        rie_condition_number: cmax / cmin.max_val(T::from_f64(1e-12)),
        cleaned_covariance: cov,
        rie_min_variance_weights: w,
        raw_min_variance_weights: mvw(&s.sample_covariance),
        rie_portfolio_variance: pv,
    }
}

#[derive(Clone, Debug)]
pub struct GarleanuPedersenDynamicAimPortfolioSpec<T: Real = f64> {
    pub covariance: Vec<Vec<T>>,
    pub current_weights: Vec<T>,
    pub slow_alpha_vector: Vec<T>,
    pub fast_alpha_vector: Vec<T>,
    pub slow_mean_reversion_phi: T,
    pub fast_mean_reversion_phi: T,
    pub risk_aversion_gamma: T,
    pub transaction_cost_kappa: T,
    pub discount_rate_r: T,
}
#[derive(Clone, Debug)]
pub struct GarleanuPedersenDynamicAimPortfolioReport<T: Real = f64> {
    pub optimal_trading_rate_a: T,
    pub slow_alpha_aim_attenuation: T,
    pub fast_alpha_aim_attenuation: T,
    pub markowitz_target_weights: Vec<T>,
    pub aim_portfolio_weights: Vec<T>,
    pub next_period_optimal_weights: Vec<T>,
    pub turnover_l1_saved_vs_markowitz: T,
    pub quadratic_transaction_cost: T,
}
pub fn solve_garleanu_pedersen_dynamic_aim_portfolio<T: Real>(
    s: &GarleanuPedersenDynamicAimPortfolioSpec<T>,
) -> GarleanuPedersenDynamicAimPortfolioReport<T> {
    let rho = T::one() / (T::one() + s.discount_rate_r.max_val(T::zero()));
    let g = s.risk_aversion_gamma.max_val(T::from_f64(1e-8));
    let k = s.transaction_cost_kappa.max_val(T::from_f64(1e-8));
    let b = g + k * (T::one() - rho);
    let a = ((b * b + T::from_f64(4.0) * g * k * rho).sqrt() - b)
        / (T::from_f64(2.0) * k * rho).max_val(T::from_f64(1e-12));
    let ps = T::one() / (T::one() + s.slow_mean_reversion_phi * k * a / g);
    let pf = T::one() / (T::one() + s.fast_mean_reversion_phi * k * a / g);
    let inv = invm(&s.covariance);
    let n = s.current_weights.len();
    let mut mw = vec![T::zero(); n];
    let mut aim = vec![T::zero(); n];
    for i in 0..n {
        for j in 0..n {
            mw[i] += inv[i][j] * (s.slow_alpha_vector[j] + s.fast_alpha_vector[j]) / g;
            aim[i] += inv[i][j] * (ps * s.slow_alpha_vector[j] + pf * s.fast_alpha_vector[j]) / g;
        }
    }
    let next = s
        .current_weights
        .iter()
        .zip(&aim)
        .map(|(x, y)| (T::one() - a) * *x + a * *y)
        .collect::<Vec<_>>();
    let markowitz_turnover = mw
        .iter()
        .zip(&s.current_weights)
        .map(|(x, y)| (*x - *y).abs())
        .fold(T::zero(), |z, x| z + x);
    let actual_turnover = next
        .iter()
        .zip(&s.current_weights)
        .map(|(x, y)| (*x - *y).abs())
        .fold(T::zero(), |z, x| z + x);
    let saved = (markowitz_turnover - actual_turnover).max_val(T::zero());
    let tc = next
        .iter()
        .zip(&s.current_weights)
        .map(|(x, y)| T::from_f64(0.5) * k * (*x - *y) * (*x - *y))
        .fold(T::zero(), |z, x| z + x);
    GarleanuPedersenDynamicAimPortfolioReport {
        optimal_trading_rate_a: a,
        slow_alpha_aim_attenuation: ps,
        fast_alpha_aim_attenuation: pf,
        markowitz_target_weights: mw,
        aim_portfolio_weights: aim,
        next_period_optimal_weights: next,
        turnover_l1_saved_vs_markowitz: saved,
        quadratic_transaction_cost: tc,
    }
}

#[derive(Clone, Debug)]
pub struct JohansenTsayVecmCointegrationSpec<T: Real = f64> {
    pub price_series: Vec<[T; 2]>,
}
#[derive(Clone, Debug)]
pub struct JohansenTsayVecmCointegrationReport<T: Real = f64> {
    pub canonical_eigenvalues: [T; 2],
    pub johansen_trace_stat_r0: T,
    pub johansen_trace_stat_r1: T,
    pub johansen_max_eigen_stat_r0: T,
    pub cointegrating_vector_beta: [T; 2],
    pub adjustment_speed_alpha: [T; 2],
    pub ou_mean_reversion_kappa: T,
    pub ou_half_life_periods: T,
    pub terminal_spread_z_score: T,
}
pub fn solve_johansen_tsay_vecm_cointegration_basket<T: Real>(
    s: &JohansenTsayVecmCointegrationSpec<T>,
) -> JohansenTsayVecmCointegrationReport<T> {
    let n = s.price_series.len().max(2);
    let teff = s.price_series.len().saturating_sub(1).max(1);
    let tf = T::from_f64(teff as f64);
    let eps = T::from_f64(1e-12);
    let mut mean_d = [T::zero(); 2];
    let mut mean_l = [T::zero(); 2];
    for t in 1..s.price_series.len() {
        mean_d[0] += s.price_series[t][0] - s.price_series[t - 1][0];
        mean_d[1] += s.price_series[t][1] - s.price_series[t - 1][1];
        mean_l[0] += s.price_series[t - 1][0];
        mean_l[1] += s.price_series[t - 1][1];
    }
    mean_d[0] /= tf;
    mean_d[1] /= tf;
    mean_l[0] /= tf;
    mean_l[1] /= tf;
    let mut s00 = [[T::zero(); 2]; 2];
    let mut s01 = [[T::zero(); 2]; 2];
    let mut s11 = [[T::zero(); 2]; 2];
    for t in 1..s.price_series.len() {
        let r0 = [
            s.price_series[t][0] - s.price_series[t - 1][0] - mean_d[0],
            s.price_series[t][1] - s.price_series[t - 1][1] - mean_d[1],
        ];
        let r1 = [
            s.price_series[t - 1][0] - mean_l[0],
            s.price_series[t - 1][1] - mean_l[1],
        ];
        for i in 0..2 {
            for j in 0..2 {
                s00[i][j] += r0[i] * r0[j] / tf;
                s01[i][j] += r0[i] * r1[j] / tf;
                s11[i][j] += r1[i] * r1[j] / tf;
            }
        }
    }
    s00[0][0] += eps;
    s00[1][1] += eps;
    s11[0][0] += eps;
    s11[1][1] += eps;
    let inv2 = |a: [[T; 2]; 2]| -> [[T; 2]; 2] {
        let mut d = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        if d.val().abs() < 1e-14 {
            d = T::from_f64(if d.val().is_sign_negative() {
                -1e-14
            } else {
                1e-14
            });
        }
        let di = T::one() / d;
        [[a[1][1] * di, -a[0][1] * di], [-a[1][0] * di, a[0][0] * di]]
    };
    let l00 = s11[0][0].sqrt();
    let l10 = s11[1][0] / l00.max_val(eps);
    let l11 = (s11[1][1] - l10 * l10).max_val(eps).sqrt();
    let linv = [
        [T::one() / l00.max_val(eps), T::zero()],
        [
            -l10 / (l00.max_val(eps) * l11.max_val(eps)),
            T::one() / l11.max_val(eps),
        ],
    ];
    let linvt = [[linv[0][0], linv[1][0]], [linv[0][1], linv[1][1]]];
    let mt = |a: [[T; 2]; 2], b: [[T; 2]; 2]| -> [[T; 2]; 2] {
        [
            [
                a[0][0] * b[0][0] + a[0][1] * b[1][0],
                a[0][0] * b[0][1] + a[0][1] * b[1][1],
            ],
            [
                a[1][0] * b[0][0] + a[1][1] * b[1][0],
                a[1][0] * b[0][1] + a[1][1] * b[1][1],
            ],
        ]
    };
    let si = inv2(s00);
    let s10 = [[s01[0][0], s01[1][0]], [s01[0][1], s01[1][1]]];
    let mtx = mt(mt(mt(linv, s10), si), mt(s01, linvt));
    let tr = mtx[0][0] + mtx[1][1];
    let disc = ((tr - T::from_f64(2.0) * mtx[0][0]).powf(2.0)
        + T::from_f64(4.0) * mtx[0][1] * mtx[0][1]
        + T::from_f64(1e-24))
    .sqrt();
    let mut l1 = ((tr + disc) / T::from_f64(2.0))
        .max_val(T::zero())
        .min_val(T::one() - eps);
    let mut l2 = ((tr - disc) / T::from_f64(2.0))
        .max_val(T::zero())
        .min_val(T::one() - eps);
    if l2.val() > l1.val() {
        std::mem::swap(&mut l1, &mut l2);
    }
    let lam = l1;
    let mut v = [mtx[0][1], lam - mtx[0][0]];
    if v[0].abs().val() + v[1].abs().val() < 1e-14 {
        v = [lam - mtx[1][1], mtx[1][0]];
    }
    let vn = (v[0] * v[0] + v[1] * v[1]).sqrt().max_val(eps);
    v[0] /= vn;
    v[1] /= vn;
    let tilde = [
        linvt[0][0] * v[0] + linvt[0][1] * v[1],
        linvt[1][0] * v[0] + linvt[1][1] * v[1],
    ];
    let t0 = if tilde[0].abs().val() < 1e-14 {
        T::from_f64(if tilde[0].val().is_sign_negative() {
            -1e-14
        } else {
            1e-14
        })
    } else {
        tilde[0]
    };
    let beta = [T::one(), tilde[1] / t0];
    let s11b = [
        s11[0][0] * beta[0] + s11[0][1] * beta[1],
        s11[1][0] * beta[0] + s11[1][1] * beta[1],
    ];
    let den = (beta[0] * s11b[0] + beta[1] * s11b[1]).max_val(eps);
    let s01b = [
        s01[0][0] * beta[0] + s01[0][1] * beta[1],
        s01[1][0] * beta[0] + s01[1][1] * beta[1],
    ];
    let alpha = [s01b[0] / den, s01b[1] / den];
    let tr0 = -tf * ((T::one() - l1).ln() + (T::one() - l2).ln());
    let tr1 = -tf * (T::one() - l2).ln();
    let max0 = -tf * (T::one() - l1).ln();
    let mut zbar = T::zero();
    for y in &s.price_series {
        zbar += y[0] + beta[1] * y[1];
    }
    zbar /= T::from_f64(n as f64);
    let mut den_ar = T::zero();
    let mut num_ar = T::zero();
    let mut vz = T::zero();
    for t in 1..s.price_series.len() {
        let z0 = s.price_series[t - 1][0] + beta[1] * s.price_series[t - 1][1] - zbar;
        let z1 = s.price_series[t][0] + beta[1] * s.price_series[t][1] - zbar;
        num_ar += z1 * z0;
        den_ar += z0 * z0;
        vz += z1 * z1;
    }
    let phi = (num_ar / den_ar.max_val(eps))
        .max_val(T::from_f64(0.01))
        .min_val(T::from_f64(0.9999));
    let kappa = -phi.ln();
    let half = T::from_f64(std::f64::consts::LN_2) / kappa.max_val(eps);
    let sigma_z = (vz / tf).sqrt();
    let last = s.price_series.last().copied().unwrap_or([T::zero(); 2]);
    let zscore = (last[0] + beta[1] * last[1] - zbar) / sigma_z.max_val(T::from_f64(1e-12));
    JohansenTsayVecmCointegrationReport {
        canonical_eigenvalues: [l1, l2],
        johansen_trace_stat_r0: tr0,
        johansen_trace_stat_r1: tr1,
        johansen_max_eigen_stat_r0: max0,
        cointegrating_vector_beta: beta,
        adjustment_speed_alpha: alpha,
        ou_mean_reversion_kappa: kappa,
        ou_half_life_periods: half,
        terminal_spread_z_score: zscore,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct LopezDePradoCausalFactorSpec<T: Real = f64> {
    pub true_causal_beta_xy: T,
    pub confounder_to_factor_axz: T,
    pub confounder_to_return_byz: T,
    pub factor_to_collider_ccx: T,
    pub return_to_collider_ccy: T,
    pub confounder_var: T,
    pub factor_idiosyncratic_var: T,
    pub return_idiosyncratic_var: T,
    pub collider_idiosyncratic_var: T,
}
#[derive(Clone, Copy, Debug)]
pub struct LopezDePradoCausalFactorReport<T: Real = f64> {
    pub unadjusted_observational_beta: T,
    pub backdoor_adjusted_beta: T,
    pub kitchen_sink_collider_conditioned_beta: T,
    pub omitted_confounder_bias: T,
    pub collider_stratification_bias: T,
    pub backdoor_recovery_abs_error: T,
}
pub fn evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias<T: Real>(
    s: &LopezDePradoCausalFactorSpec<T>,
) -> LopezDePradoCausalFactorReport<T> {
    let vz = s.confounder_var;
    let cxz = s.confounder_to_factor_axz * vz;
    let vx =
        s.confounder_to_factor_axz * s.confounder_to_factor_axz * vz + s.factor_idiosyncratic_var;
    let cyz = s.true_causal_beta_xy * cxz + s.confounder_to_return_byz * vz;
    let cxy = s.true_causal_beta_xy * vx + s.confounder_to_return_byz * cxz;
    let vy = s.true_causal_beta_xy * s.true_causal_beta_xy * vx
        + s.confounder_to_return_byz * s.confounder_to_return_byz * vz
        + T::from_f64(2.0) * s.true_causal_beta_xy * s.confounder_to_return_byz * cxz
        + s.return_idiosyncratic_var;
    let naive = cxy / vx.max_val(T::from_f64(1e-14));
    let cxy_z = cxy - cxz * cyz / vz.max_val(T::from_f64(1e-14));
    let vx_z = vx - cxz * cxz / vz.max_val(T::from_f64(1e-14));
    let vy_z = vy - cyz * cyz / vz.max_val(T::from_f64(1e-14));
    let back = cxy_z / vx_z.max_val(T::from_f64(1e-14));
    let xcz = s.factor_to_collider_ccx * vx_z + s.return_to_collider_ccy * cxy_z;
    let ycz = s.factor_to_collider_ccx * cxy_z + s.return_to_collider_ccy * vy_z;
    let vcz = s.factor_to_collider_ccx * s.factor_to_collider_ccx * vx_z
        + s.return_to_collider_ccy * s.return_to_collider_ccy * vy_z
        + T::from_f64(2.0) * s.factor_to_collider_ccx * s.return_to_collider_ccy * cxy_z
        + s.collider_idiosyncratic_var;
    let cxy_zc = cxy_z - xcz * ycz / vcz.max_val(T::from_f64(1e-14));
    let vx_zc = vx_z - xcz * xcz / vcz.max_val(T::from_f64(1e-14));
    let collider = cxy_zc / vx_zc.max_val(T::from_f64(1e-14));
    LopezDePradoCausalFactorReport {
        unadjusted_observational_beta: naive,
        backdoor_adjusted_beta: back,
        kitchen_sink_collider_conditioned_beta: collider,
        omitted_confounder_bias: naive - s.true_causal_beta_xy,
        collider_stratification_bias: collider - s.true_causal_beta_xy,
        backdoor_recovery_abs_error: (back - s.true_causal_beta_xy).abs(),
    }
}

#[derive(Clone, Copy, Debug)]
pub struct MadhavanAbnerHillEtfApArbitrageSpec<T: Real = f64> {
    pub etf_market_price: T,
    pub reported_inav: T,
    pub futures_beta_return: T,
    pub stale_underlying_weight: T,
    pub ap_creation_fee_bps: T,
    pub basket_bid_ask_spread_bps: T,
    pub creation_unit_volatility_bps: T,
    pub ap_risk_aversion: T,
    pub embedded_unrealized_gain_fraction: T,
    pub annual_redemption_turnover: T,
    pub capital_gains_tax_rate: T,
}
#[derive(Clone, Copy, Debug)]
pub struct MadhavanAbnerHillEtfApArbitrageReport<T: Real = f64> {
    pub synchronous_economic_nav: T,
    pub illusionary_inav_premium_bps: T,
    pub true_arbitrage_basis_bps: T,
    pub ap_no_arbitrage_half_band_bps: T,
    pub net_ap_creation_redemption_profit_bps: T,
    pub optimal_ap_creation_intensity: T,
    pub annual_in_kind_heartbeat_tax_alpha_bps: T,
    pub etf_vs_mutual_fund_after_tax_advantage_bps: T,
}
pub fn evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha<T: Real>(
    s: &MadhavanAbnerHillEtfApArbitrageSpec<T>,
) -> MadhavanAbnerHillEtfApArbitrageReport<T> {
    let sync = s.reported_inav * (T::one() - s.stale_underlying_weight)
        + s.reported_inav * (T::one() + s.futures_beta_return) * s.stale_underlying_weight;
    let illusion =
        (s.etf_market_price / s.reported_inav.max_val(T::from_f64(1e-10)) - T::one()) * 10000.0;
    let basis = (s.etf_market_price / sync.max_val(T::from_f64(1e-10)) - T::one()) * 10000.0;
    let band = s.ap_creation_fee_bps
        + s.basket_bid_ask_spread_bps / 2.0
        + s.creation_unit_volatility_bps * s.ap_risk_aversion;
    let profit = (basis.abs() - band).max_val(T::zero());
    let intensity = (profit / band.max_val(T::from_f64(1e-8))).tanh();
    let tax = s.embedded_unrealized_gain_fraction
        * s.annual_redemption_turnover
        * s.capital_gains_tax_rate
        * 10000.0;
    MadhavanAbnerHillEtfApArbitrageReport {
        synchronous_economic_nav: sync,
        illusionary_inav_premium_bps: illusion,
        true_arbitrage_basis_bps: basis,
        ap_no_arbitrage_half_band_bps: band,
        net_ap_creation_redemption_profit_bps: profit,
        optimal_ap_creation_intensity: intensity,
        annual_in_kind_heartbeat_tax_alpha_bps: tax.max_val(T::zero()),
        etf_vs_mutual_fund_after_tax_advantage_bps: tax
            - s.annual_redemption_turnover * s.ap_creation_fee_bps,
    }
}

#[derive(Clone, Debug)]
pub struct GrinoldQianChristophersonAlphaDecaySpec<T: Real = f64> {
    pub fast_signal_ic0: T,
    pub fast_signal_half_life_months: T,
    pub slow_signal_ic0: T,
    pub slow_signal_half_life_months: T,
    pub signal_correlation: T,
    pub rebalance_horizon_months: T,
    pub transaction_cost_penalty: T,
    pub factor_active_exposures: Vec<T>,
    pub factor_benchmark_returns: Vec<T>,
    pub normal_factor_exposures: Vec<T>,
}
#[derive(Clone, Debug)]
pub struct GrinoldQianChristophersonAlphaDecayReport<T: Real = f64> {
    pub fast_decayed_ic: T,
    pub slow_decayed_ic: T,
    pub optimal_fast_signal_weight: T,
    pub optimal_slow_signal_weight: T,
    pub composite_blended_ic: T,
    pub turnover_adjusted_information_ratio: T,
    pub christopherson_policy_return: T,
    pub christopherson_factor_timing_return: T,
    pub total_factor_contribution: T,
}
pub fn evaluate_grinold_qian_christopherson_ic_decay_and_timing<T: Real>(
    s: &GrinoldQianChristophersonAlphaDecaySpec<T>,
) -> GrinoldQianChristophersonAlphaDecayReport<T> {
    let fd = s.fast_signal_ic0
        * (-(T::from_f64(2.0).ln()) * s.rebalance_horizon_months
            / s.fast_signal_half_life_months.max_val(T::from_f64(1e-8)))
        .exp();
    let sd = s.slow_signal_ic0
        * (-(T::from_f64(2.0).ln()) * s.rebalance_horizon_months
            / s.slow_signal_half_life_months.max_val(T::from_f64(1e-8)))
        .exp();
    let kf = T::from_f64(2.0).ln() / s.fast_signal_half_life_months.max_val(T::from_f64(1e-8));
    let ks = T::from_f64(2.0).ln() / s.slow_signal_half_life_months.max_val(T::from_f64(1e-8));
    let rho = s.signal_correlation;
    let af = T::one() + s.transaction_cost_penalty * ks;
    let as_ = T::one() + s.transaction_cost_penalty * kf;
    let den = (af * as_ - rho * rho).max_val(T::from_f64(1e-8));
    let twf = (fd * af - rho * sd) / den;
    let tws = (sd * as_ - rho * fd) / den;
    let twsum = (twf + tws).max_val(T::from_f64(1e-8));
    let wf = twf / twsum;
    let ws = T::one() - wf;
    let norm = (wf * wf + ws * ws + T::from_f64(2.0) * rho * wf * ws)
        .max_val(T::from_f64(1e-12))
        .sqrt();
    let blend = (wf * fd + ws * sd) / norm;
    let ir = blend / (T::one() + s.transaction_cost_penalty.abs());
    let mut policy = T::zero();
    let mut timing = T::zero();
    for i in 0..s
        .factor_active_exposures
        .len()
        .min(s.factor_benchmark_returns.len())
        .min(s.normal_factor_exposures.len())
    {
        policy += s.normal_factor_exposures[i] * s.factor_benchmark_returns[i];
        timing += (s.factor_active_exposures[i] - s.normal_factor_exposures[i])
            * s.factor_benchmark_returns[i];
    }
    let total = policy + timing;
    GrinoldQianChristophersonAlphaDecayReport {
        fast_decayed_ic: fd,
        slow_decayed_ic: sd,
        optimal_fast_signal_weight: wf,
        optimal_slow_signal_weight: ws,
        composite_blended_ic: blend,
        turnover_adjusted_information_ratio: ir,
        christopherson_policy_return: policy,
        christopherson_factor_timing_return: timing,
        total_factor_contribution: total,
    }
}

#[cfg(test)]
mod phase23_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn nco_weights_sum() {
        let s = LopezDePradoNcoDenoisedSpec {
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            expected_returns: vec![0.05, 0.06],
            sample_ratio_q: 2.0,
            denoising_alpha: 0.5,
            num_clusters: 2,
        };
        let r = solve_lopez_de_prado_denoised_nco_allocation(&s);
        assert!((r.nco_optimal_weights.iter().sum::<f64>() - 1.0).abs() < 1e-8);
        assert!((r.nco_optimal_weights[0] - r.nco_optimal_weights[1]).abs() > 1e-8);
        assert!(r.denoised_condition_number <= r.raw_condition_number + 1e-8);
    }

    #[test]
    fn dual_fd_parity() {
        fn mk<T: Real>(x: T) -> LopezDePradoNcoDenoisedSpec<T> {
            LopezDePradoNcoDenoisedSpec {
                covariance: vec![
                    vec![T::from_f64(0.04), T::from_f64(0.01)],
                    vec![T::from_f64(0.01), T::from_f64(0.09)],
                ],
                expected_returns: vec![x, T::from_f64(0.06)],
                sample_ratio_q: T::from_f64(2.0),
                denoising_alpha: T::from_f64(0.5),
                num_clusters: 2,
            }
        }
        let h = 1e-5;
        let fd = (solve_lopez_de_prado_denoised_nco_allocation(&mk(0.05 + h))
            .nco_portfolio_volatility
            .val()
            - solve_lopez_de_prado_denoised_nco_allocation(&mk(0.05 - h))
                .nco_portfolio_volatility
                .val())
            / (2.0 * h);
        let ad = solve_lopez_de_prado_denoised_nco_allocation(&mk(Dual::<1>::variable(0.05, 0)))
            .nco_portfolio_volatility;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase14_tests {
    use super::*;
    #[test]
    fn phase14_hrp_diagonal_inverse_variance() {
        let c = vec![vec![0.04, 0.0], vec![0.0, 0.09]];
        let w = hierarchical_risk_parity_weights(&c);
        assert!((w[0] - 0.09 / 0.13).abs() < 1e-12);
        assert!((w[1] - 0.04 / 0.13).abs() < 1e-12);
    }
}

#[cfg(test)]
mod phase13_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn robust_markowitz_dual_matches_fd_in_mu_and_kappa() {
        let eval_f64 = |mu0: f64, kappa: f64| {
            robust_markowitz_weights(&RobustMarkowitzSpec {
                covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
                return_uncertainty_cov: vec![vec![0.02, 0.005], vec![0.005, 0.03]],
                expected_returns: vec![mu0, 0.04],
                risk_aversion: 2.5,
                uncertainty_radius: kappa,
                lower_bounds: vec![0.0, 0.0],
                upper_bounds: vec![1.0, 1.0],
            })[0]
        };
        let w_ad = robust_markowitz_weights(&RobustMarkowitzSpec {
            covariance: vec![
                vec![Dual::<2>::from_value(0.04), Dual::<2>::from_value(0.01)],
                vec![Dual::<2>::from_value(0.01), Dual::<2>::from_value(0.09)],
            ],
            return_uncertainty_cov: vec![
                vec![Dual::<2>::from_value(0.02), Dual::<2>::from_value(0.005)],
                vec![Dual::<2>::from_value(0.005), Dual::<2>::from_value(0.03)],
            ],
            expected_returns: vec![Dual::<2>::variable(0.07, 0), Dual::<2>::from_value(0.04)],
            risk_aversion: Dual::<2>::from_value(2.5),
            uncertainty_radius: Dual::<2>::variable(0.15, 1),
            lower_bounds: vec![Dual::<2>::from_value(0.0), Dual::<2>::from_value(0.0)],
            upper_bounds: vec![Dual::<2>::from_value(1.0), Dual::<2>::from_value(1.0)],
        })[0];
        let h = 1e-6;
        let fd_mu0 = (eval_f64(0.07 + h, 0.15) - eval_f64(0.07 - h, 0.15)) / (2.0 * h);
        let fd_kap = (eval_f64(0.07, 0.15 + h) - eval_f64(0.07, 0.15 - h)) / (2.0 * h);
        assert!((w_ad.d[0] - fd_mu0).abs() < 1e-6);
        assert!((w_ad.d[1] - fd_kap).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase15_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn gerber_rmt_nco() {
        let r = vec![
            vec![0.01, 0.02],
            vec![-0.02, -0.01],
            vec![0.03, 0.02],
            vec![-0.01, -0.02],
        ];
        let g = gerber_covariance_matrix(&r, 0.5, 1e-3).unwrap();
        assert!(g[0][0] > 0.0 && (g[0][1] - g[1][0]).abs() < 1e-12);
        let c = vec![
            vec![0.04, 0.01, 0., 0.],
            vec![0.01, 0.09, 0., 0.],
            vec![0., 0., 0.16, 0.02],
            vec![0., 0., 0.02, 0.25],
        ];
        let w = nested_clustered_optimization_weights(&NestedClusteredOptimizationSpec {
            covariance: c.clone(),
            expected_returns: None,
            risk_aversion: 1.0,
            num_clusters: 2,
        });
        let mv = markowitz_minimum_variance_weights(&c);
        assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-10);
        for i in 0..4 {
            assert!((w[i] - mv[i]).abs() < 1e-10);
        }
        let h = 1e-5;
        let eval_rmt = |x01: f64| {
            let mut c2 = c.clone();
            c2[0][1] = x01;
            c2[1][0] = x01;
            marchenko_pastur_denoise_covariance(&c2, 12).unwrap()[0][1]
        };
        let fd_rmt = (eval_rmt(0.01 + h) - eval_rmt(0.01 - h)) / (2.0 * h);
        type D = Dual<1>;
        let mut cd: Vec<Vec<D>> = c
            .iter()
            .map(|row| row.iter().map(|&v| D::from_value(v)).collect())
            .collect();
        cd[0][1] = D::variable(0.01, 0);
        cd[1][0] = D::variable(0.01, 0);
        let rmt_d = marchenko_pastur_denoise_covariance(&cd, 12).unwrap();
        assert!((rmt_d[0][1].d[0] - fd_rmt).abs() < 1e-5);
        let eval_nco0 = |x01: f64| {
            let mut c2 = c.clone();
            c2[0][1] = x01;
            c2[1][0] = x01;
            nested_clustered_optimization_weights(&NestedClusteredOptimizationSpec {
                covariance: c2,
                expected_returns: None,
                risk_aversion: 1.0,
                num_clusters: 2,
            })[0]
        };
        let fd_nco = (eval_nco0(0.01 + h) - eval_nco0(0.01 - h)) / (2.0 * h);
        let nco_d = nested_clustered_optimization_weights(&NestedClusteredOptimizationSpec {
            covariance: cd,
            expected_returns: None,
            risk_aversion: D::one(),
            num_clusters: 2,
        });
        assert!((nco_d[0].d[0] - fd_nco).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase17_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn idzorek_boundary_limits_and_kelly_ad() {
        let cov = vec![vec![0.04, 0.01], vec![0.01, 0.09]];
        let w_mkt = vec![0.6, 0.4];
        let mut s = IdzorekBlackLittermanSpec {
            covariance: cov.clone(),
            market_weights: w_mkt.clone(),
            risk_aversion: 2.5,
            tau: 0.05,
            view_matrix: vec![vec![1.0, -1.0]],
            view_returns: vec![0.03],
            view_confidences: vec![0.0],
        };
        let p0 = idzorek_black_litterman_posterior(&s);
        assert!((p0.weights[0] - w_mkt[0]).abs() < 1e-12);
        assert!((p0.weights[1] - w_mkt[1]).abs() < 1e-12);

        s.view_confidences = vec![1.0];
        let p1 = idzorek_black_litterman_posterior(&s);
        let view_realized = p1.expected_returns[0] - p1.expected_returns[1];
        assert!((view_realized - 0.03).abs() < 1e-12);

        let h = 1e-5;
        let eval_c = |c_val: f64| {
            let mut sp = s.clone();
            sp.view_confidences = vec![c_val];
            idzorek_black_litterman_posterior(&sp).weights[0]
        };
        let fd = (eval_c(0.65 + h) - eval_c(0.65 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = IdzorekBlackLittermanSpec {
            covariance: vec![
                vec![D::from_value(0.04), D::from_value(0.01)],
                vec![D::from_value(0.01), D::from_value(0.09)],
            ],
            market_weights: vec![D::from_value(0.6), D::from_value(0.4)],
            risk_aversion: D::from_value(2.5),
            tau: D::from_value(0.05),
            view_matrix: vec![vec![D::one(), -D::one()]],
            view_returns: vec![D::from_value(0.03)],
            view_confidences: vec![D::variable(0.65, 0)],
        };
        let ad = idzorek_black_litterman_posterior(&sd).weights[0];
        assert!((ad.d[0] - fd).abs() < 1e-6);

        let kelly = kelly_fractional_allocation(&KellyAllocationSpec {
            expected_returns: vec![0.08, 0.11],
            covariance: cov,
            risk_free_rate: 0.03,
            max_kelly_fraction: 0.5,
            max_drawdown_limit: 0.20,
            drawdown_tail_prob: 0.10,
        });
        let kelly_tighter = kelly_fractional_allocation(&KellyAllocationSpec {
            expected_returns: vec![0.08, 0.11],
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            risk_free_rate: 0.03,
            max_kelly_fraction: 0.5,
            max_drawdown_limit: 0.10,
            drawdown_tail_prob: 0.10,
        });
        assert!(kelly.effective_kelly_fraction > 0.0 && kelly.effective_kelly_fraction < 0.5);
        assert!(kelly_tighter.effective_kelly_fraction < kelly.effective_kelly_fraction);
        assert!(kelly.expected_log_growth_rate > 0.03);
    }
}

#[cfg(test)]
mod phase18_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn entropy_pooling_view_satisfaction_weights_and_ad() {
        let s_prior = EntropyPoolingSpec {
            scenario_returns: vec![vec![0.08, 0.02], vec![-0.04, 0.06], vec![0.02, -0.03]],
            prior_probs: vec![1.0 / 3.0; 3],
            view_matrix: vec![vec![0.08, -0.04, 0.02]],
            view_targets: vec![0.02],
            risk_aversion: 2.5,
        };
        let r0 = entropy_pooling_posterior(&s_prior);
        assert!(r0.kl_divergence.abs() < 1e-12);

        let mut s_view = s_prior.clone();
        s_view.view_targets = vec![0.04];
        let r1 = entropy_pooling_posterior(&s_view);
        assert!(r1.kl_divergence > 0.0);
        assert!((r1.posterior_mean[0] - 0.04).abs() < 1e-10);
        assert!(r1.weights[0].abs() > 0.0);

        let h = 1e-5;
        let eval_b = |b: f64| {
            let mut sp = s_view.clone();
            sp.view_targets = vec![b];
            entropy_pooling_posterior(&sp).weights[0]
        };
        let fd = (eval_b(0.04 + h) - eval_b(0.04 - h)) / (2.0 * h);
        type D = Dual<1>;
        let sd = EntropyPoolingSpec {
            scenario_returns: vec![
                vec![D::from_value(0.08), D::from_value(0.02)],
                vec![D::from_value(-0.04), D::from_value(0.06)],
                vec![D::from_value(0.02), D::from_value(-0.03)],
            ],
            prior_probs: vec![D::from_value(1.0 / 3.0); 3],
            view_matrix: vec![vec![
                D::from_value(0.08),
                D::from_value(-0.04),
                D::from_value(0.02),
            ]],
            view_targets: vec![D::variable(0.04, 0)],
            risk_aversion: D::from_value(2.5),
        };
        let ad = entropy_pooling_posterior(&sd).weights[0];
        assert!((ad.d[0] - fd).abs() < 5e-5);
    }
}

#[cfg(test)]
mod phase19_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn cotton_schur_hrp_unifies_hrp_and_minimum_variance_and_oas_shrinkage_ad() {
        // Cotton (2024) Theorem: on a 3x3 equicorrelated unit-variance matrix,
        // standard HRP (gamma = 0) gives asymmetric weights, whereas Schur-HRP at gamma = 1
        // recovers the exact symmetric Minimum-Variance weights (1/3, 1/3, 1/3)!
        let rho = 0.5;
        let cov = vec![
            vec![1.0, rho, rho],
            vec![rho, 1.0, rho],
            vec![rho, rho, 1.0],
        ];
        let hrp0 = schur_complementary_hrp(&SchurHrpSpec {
            covariance: cov.clone(),
            gamma: 0.0,
        });
        let mv1 = schur_complementary_hrp(&SchurHrpSpec {
            covariance: cov.clone(),
            gamma: 1.0,
        });
        assert!((hrp0.weights[0] - 1.0 / 3.0).abs() > 1e-2);
        for &w in &mv1.weights {
            assert!((w - 1.0 / 3.0).abs() < 1e-12);
        }
        assert!(mv1.portfolio_variance < hrp0.portfolio_variance);

        let rets = vec![
            vec![0.02, -0.01, 0.015],
            vec![-0.01, 0.03, -0.005],
            vec![0.015, 0.005, 0.02],
            vec![-0.02, -0.015, 0.01],
            vec![0.01, 0.02, -0.01],
        ];
        let shrink = oas_and_nonlinear_shrinkage_covariance(&rets);
        assert!(shrink.oas_shrinkage_intensity > 0.0 && shrink.oas_shrinkage_intensity <= 1.0);
        let tr_oas =
            shrink.oas_covariance[0][0] + shrink.oas_covariance[1][1] + shrink.oas_covariance[2][2];
        let tr_nl = shrink.nonlinear_shrunk_covariance[0][0]
            + shrink.nonlinear_shrunk_covariance[1][1]
            + shrink.nonlinear_shrunk_covariance[2][2];
        assert!((tr_oas - tr_nl).abs() < 1e-12);

        let h = 1e-5;
        let eval_g = |g: f64| {
            schur_complementary_hrp(&SchurHrpSpec {
                covariance: cov.clone(),
                gamma: g,
            })
            .weights[0]
        };
        let fd = (eval_g(0.5 + h) - eval_g(0.5 - h)) / (2.0 * h);
        type D = Dual<1>;
        let cov_d: Vec<Vec<D>> = cov
            .iter()
            .map(|r| r.iter().map(|&x| D::from_value(x)).collect())
            .collect();
        let ad = schur_complementary_hrp(&SchurHrpSpec {
            covariance: cov_d,
            gamma: D::variable(0.5, 0),
        })
        .weights[0];
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase20_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn principal_component_risk_parity_exact_factor_shares_enb_and_ad() {
        let cov = vec![
            vec![0.04, 0.025, 0.010],
            vec![0.025, 0.05, 0.015],
            vec![0.010, 0.015, 0.03],
        ];
        let spec = FactorEigenRiskParitySpec {
            covariance: cov.clone(),
            factor_risk_budgets: vec![1.0 / 3.0; 3],
        };
        let rep = solve_principal_component_risk_parity(&spec);
        assert!((rep.weights.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        // Equal factor risk budgets => each principal component contributes 1/3 of variance and PCA ENB == 3.0!
        for &s in &rep.factor_variance_shares {
            assert!((s - 1.0 / 3.0).abs() < 1e-10);
        }
        assert!((rep.effective_number_of_pca_bets - 3.0).abs() < 1e-10);
        assert!(rep.effective_number_of_pca_bets > rep.equal_weight_pca_enb);

        let h = 1e-5;
        let eval_b0 = |b0: f64| {
            solve_principal_component_risk_parity(&FactorEigenRiskParitySpec {
                covariance: cov.clone(),
                factor_risk_budgets: vec![b0, 1.0 / 3.0, 1.0 / 3.0],
            })
            .weights[0]
        };
        let fd = (eval_b0(1.0 / 3.0 + h) - eval_b0(1.0 / 3.0 - h)) / (2.0 * h);
        type D = Dual<1>;
        let cov_d: Vec<Vec<D>> = cov
            .iter()
            .map(|r| r.iter().map(|&x| D::from_value(x)).collect())
            .collect();
        let ad = solve_principal_component_risk_parity(&FactorEigenRiskParitySpec {
            covariance: cov_d,
            factor_risk_budgets: vec![
                D::variable(1.0 / 3.0, 0),
                D::from_value(1.0 / 3.0),
                D::from_value(1.0 / 3.0),
            ],
        })
        .weights[0];
        assert!((ad.d[0] - fd).abs() < 1e-6);
    }
}

#[cfg(test)]
mod phase21_tests {
    use super::*;

    #[test]
    fn paleologo_factor_neutrality_and_crowding() {
        let s = PaleologoFactorNeutralAlphaSpec {
            raw_alphas: vec![1.0, 2.0, 1.5],
            factor_loadings_b: vec![vec![1.0], vec![1.0], vec![0.5]],
            factor_covariance_omega: vec![vec![1.0]],
            sample_idio_variances: vec![0.1, 0.2, 0.15],
            idio_shrinkage_xi: 0.5,
            idio_risk_aversion: 1.0,
            factor_risk_aversion: 2.0,
            factor_crowding_scores: vec![1.0],
            crowding_penalty_lambda: 1.0,
        };
        let r = solve_paleologo_factor_neutral_portfolio(&s);
        assert!(r.residual_factor_exposures[0].val().is_finite());
        assert!(
            r.crowding_penalized_weights
                .iter()
                .all(|x| x.val().is_finite())
        );
    }
}

#[cfg(test)]
mod phase22_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn meucci_solver_finite() {
        let s = MeucciIdzorekBlackLittermanSpec {
            covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            market_weights: vec![0.6, 0.4],
            risk_aversion_lambda: 2.5,
            tau: 0.05,
            view_matrix_p: vec![vec![1.0, -1.0]],
            view_returns_q: vec![0.03],
            view_confidences: vec![0.65],
        };
        let r = solve_meucci_idzorek_black_litterman_and_plausibility(&s);
        assert!(r.mahalanobis_plausibility_squared.val().is_finite());
        let h = 1e-5;
        let mut up = s.clone();
        up.view_returns_q[0] += h;
        let mut dn = s.clone();
        dn.view_returns_q[0] -= h;
        let fd = (solve_meucci_idzorek_black_litterman_and_plausibility(&up)
            .posterior_expected_returns[0]
            - solve_meucci_idzorek_black_litterman_and_plausibility(&dn)
                .posterior_expected_returns[0])
            / (2.0 * h);
        let ad = solve_meucci_idzorek_black_litterman_and_plausibility(
            &MeucciIdzorekBlackLittermanSpec {
                covariance: vec![
                    vec![Dual::<1>::from_value(0.04), Dual::<1>::from_value(0.01)],
                    vec![Dual::<1>::from_value(0.01), Dual::<1>::from_value(0.09)],
                ],
                market_weights: vec![Dual::<1>::from_value(0.6), Dual::<1>::from_value(0.4)],
                risk_aversion_lambda: Dual::<1>::from_value(2.5),
                tau: Dual::<1>::from_value(0.05),
                view_matrix_p: vec![vec![
                    Dual::<1>::from_value(1.0),
                    Dual::<1>::from_value(-1.0),
                ]],
                view_returns_q: vec![Dual::<1>::variable(0.03, 0)],
                view_confidences: vec![Dual::<1>::from_value(0.65)],
            },
        );
        assert!((ad.posterior_expected_returns[0].d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase24_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn rie_and_gp_phase24() {
        let r = solve_potters_bouchaud_rie_nonlinear_shrinkage(&PottersBouchaudRieCovarianceSpec {
            sample_covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
            num_observations_t: 100,
        });
        assert!(
            (r.rie_cleaned_eigenvalues
                .iter()
                .copied()
                .fold(0.0, |a, x| a + x)
                - r.raw_eigenvalues.iter().copied().fold(0.0, |a, x| a + x))
            .abs()
                < 1e-11
        );
        assert!(r.rie_condition_number < r.raw_condition_number);
        fn mk<T: Real>(x: T) -> PottersBouchaudRieCovarianceSpec<T> {
            PottersBouchaudRieCovarianceSpec {
                sample_covariance: vec![
                    vec![x, T::from_f64(0.01)],
                    vec![T::from_f64(0.01), T::from_f64(0.09)],
                ],
                num_observations_t: 100,
            }
        }
        let h = 1e-5;
        let fd = (solve_potters_bouchaud_rie_nonlinear_shrinkage(&mk(0.04 + h))
            .rie_portfolio_variance
            - solve_potters_bouchaud_rie_nonlinear_shrinkage(&mk(0.04 - h)).rie_portfolio_variance)
            / (2.0 * h);
        let ad = solve_potters_bouchaud_rie_nonlinear_shrinkage(&mk(Dual::<1>::variable(0.04, 0)))
            .rie_portfolio_variance;
        assert!((ad.d[0] - fd).abs() < 1e-5);
        let g = solve_garleanu_pedersen_dynamic_aim_portfolio(
            &GarleanuPedersenDynamicAimPortfolioSpec {
                covariance: vec![
                    vec![Dual::<1>::from_value(0.04), Dual::<1>::from_value(0.01)],
                    vec![Dual::<1>::from_value(0.01), Dual::<1>::from_value(0.09)],
                ],
                current_weights: vec![Dual::<1>::from_value(0.5), Dual::<1>::from_value(0.5)],
                slow_alpha_vector: vec![Dual::<1>::from_value(0.03), Dual::<1>::from_value(0.02)],
                fast_alpha_vector: vec![Dual::<1>::from_value(0.01), Dual::<1>::from_value(0.02)],
                slow_mean_reversion_phi: Dual::<1>::from_value(0.1),
                fast_mean_reversion_phi: Dual::<1>::from_value(0.5),
                risk_aversion_gamma: Dual::<1>::from_value(2.0),
                transaction_cost_kappa: Dual::<1>::variable(0.1, 0),
                discount_rate_r: Dual::<1>::from_value(0.02),
            },
        );
        assert!(g.optimal_trading_rate_a.val() > 0.0 && g.optimal_trading_rate_a.val() < 1.0);
        assert!(g.fast_alpha_aim_attenuation < g.slow_alpha_aim_attenuation);
        assert!(g.turnover_l1_saved_vs_markowitz >= Dual::<1>::zero());
        let fd_a = (solve_garleanu_pedersen_dynamic_aim_portfolio(
            &GarleanuPedersenDynamicAimPortfolioSpec {
                covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
                current_weights: vec![0.5, 0.5],
                slow_alpha_vector: vec![0.03, 0.02],
                fast_alpha_vector: vec![0.01, 0.02],
                slow_mean_reversion_phi: 0.1,
                fast_mean_reversion_phi: 0.5,
                risk_aversion_gamma: 2.0,
                transaction_cost_kappa: 0.1 + 1e-5,
                discount_rate_r: 0.02,
            },
        )
        .optimal_trading_rate_a
            - solve_garleanu_pedersen_dynamic_aim_portfolio(
                &GarleanuPedersenDynamicAimPortfolioSpec {
                    covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
                    current_weights: vec![0.5, 0.5],
                    slow_alpha_vector: vec![0.03, 0.02],
                    fast_alpha_vector: vec![0.01, 0.02],
                    slow_mean_reversion_phi: 0.1,
                    fast_mean_reversion_phi: 0.5,
                    risk_aversion_gamma: 2.0,
                    transaction_cost_kappa: 0.1 - 1e-5,
                    discount_rate_r: 0.02,
                },
            )
            .optimal_trading_rate_a)
            / (2e-5);
        assert!((g.optimal_trading_rate_a.d[0] - fd_a).abs() < 1e-5);
        let gp_f64 = |k: f64| {
            solve_garleanu_pedersen_dynamic_aim_portfolio(
                &GarleanuPedersenDynamicAimPortfolioSpec {
                    covariance: vec![vec![0.04, 0.01], vec![0.01, 0.09]],
                    current_weights: vec![0.5, 0.5],
                    slow_alpha_vector: vec![0.03, 0.02],
                    fast_alpha_vector: vec![0.01, 0.02],
                    slow_mean_reversion_phi: 0.1,
                    fast_mean_reversion_phi: 0.5,
                    risk_aversion_gamma: 2.0,
                    transaction_cost_kappa: k,
                    discount_rate_r: 0.02,
                },
            )
        };
        let fd_next = (gp_f64(0.1 + 1e-5).next_period_optimal_weights[0]
            - gp_f64(0.1 - 1e-5).next_period_optimal_weights[0])
            / (2e-5);
        let ad_next = solve_garleanu_pedersen_dynamic_aim_portfolio(
            &GarleanuPedersenDynamicAimPortfolioSpec {
                covariance: vec![
                    vec![Dual::<1>::from_value(0.04), Dual::<1>::from_value(0.01)],
                    vec![Dual::<1>::from_value(0.01), Dual::<1>::from_value(0.09)],
                ],
                current_weights: vec![Dual::<1>::from_value(0.5), Dual::<1>::from_value(0.5)],
                slow_alpha_vector: vec![Dual::<1>::from_value(0.03), Dual::<1>::from_value(0.02)],
                fast_alpha_vector: vec![Dual::<1>::from_value(0.01), Dual::<1>::from_value(0.02)],
                slow_mean_reversion_phi: Dual::<1>::from_value(0.1),
                fast_mean_reversion_phi: Dual::<1>::from_value(0.5),
                risk_aversion_gamma: Dual::<1>::from_value(2.0),
                transaction_cost_kappa: Dual::<1>::variable(0.1, 0),
                discount_rate_r: Dual::<1>::from_value(0.02),
            },
        )
        .next_period_optimal_weights[0];
        assert!((ad_next.d[0] - fd_next).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase25_tests {
    use super::*;
    use ql_ad::Dual;
    fn causal<T: Real>(x: T) -> LopezDePradoCausalFactorSpec<T> {
        LopezDePradoCausalFactorSpec {
            true_causal_beta_xy: x,
            confounder_to_factor_axz: T::from_f64(0.8),
            confounder_to_return_byz: T::from_f64(0.6),
            factor_to_collider_ccx: T::from_f64(0.4),
            return_to_collider_ccy: T::from_f64(0.5),
            confounder_var: T::one(),
            factor_idiosyncratic_var: T::from_f64(0.5),
            return_idiosyncratic_var: T::from_f64(0.7),
            collider_idiosyncratic_var: T::from_f64(0.4),
        }
    }
    #[test]
    fn causal_backdoor_dual_fd() {
        let r0 = evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias(&causal(0.0));
        assert!(r0.backdoor_recovery_abs_error < 1e-12);
        assert!(r0.collider_stratification_bias < 0.0);
        let h = 1e-5;
        let fd = (evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias(&causal(0.5 + h))
            .unadjusted_observational_beta
            .val()
            - evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias(&causal(0.5 - h))
                .unadjusted_observational_beta
                .val())
            / (2.0 * h);
        let ad_naive = evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias(&causal(
            Dual::<1>::variable(0.5, 0),
        ))
        .unadjusted_observational_beta;
        assert!((ad_naive.d[0] - fd).abs() < 1e-5);
        let ad = evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias(&causal(
            Dual::<1>::variable(0.5, 0),
        ))
        .kitchen_sink_collider_conditioned_beta;
        let fd_collider =
            (evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias(&causal(0.5 + h))
                .kitchen_sink_collider_conditioned_beta
                .val()
                - evaluate_lopez_de_prado_causal_backdoor_vs_collider_bias(&causal(0.5 - h))
                    .kitchen_sink_collider_conditioned_beta
                    .val())
                / (2.0 * h);
        assert!((ad.d[0] - fd_collider).abs() < 1e-5);
    }
    fn vecm<T: Real>(x: T) -> JohansenTsayVecmCointegrationSpec<T> {
        JohansenTsayVecmCointegrationSpec {
            price_series: vec![
                [T::from_f64(100.0), T::from_f64(100.5) + x],
                [T::from_f64(101.0), T::from_f64(100.8) + x],
                [T::from_f64(99.5), T::from_f64(100.2) + x],
                [T::from_f64(100.8), T::from_f64(100.4) + x],
                [T::from_f64(100.2), T::from_f64(100.6) + x],
                [T::from_f64(101.2), T::from_f64(100.9) + x],
            ],
        }
    }
    #[test]
    fn vecm_dual_fd_and_trace_partition() {
        let r = solve_johansen_tsay_vecm_cointegration_basket(&vecm(0.0));
        assert!(
            (r.johansen_trace_stat_r0 - (r.johansen_max_eigen_stat_r0 + r.johansen_trace_stat_r1))
                .abs()
                < 1e-12
        );
        assert!(r.ou_mean_reversion_kappa > 0.0);
        assert!(r.ou_half_life_periods > 0.0);
        let h = 1e-5;
        let fd = (solve_johansen_tsay_vecm_cointegration_basket(&vecm(h))
            .terminal_spread_z_score
            .val()
            - solve_johansen_tsay_vecm_cointegration_basket(&vecm(-h))
                .terminal_spread_z_score
                .val())
            / (2.0 * h);
        let ad = solve_johansen_tsay_vecm_cointegration_basket(&vecm(Dual::<1>::variable(0.0, 0)))
            .terminal_spread_z_score;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase26_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(ic: T) -> GrinoldKahnClarkeFundamentalLawSpec<T> {
        GrinoldKahnClarkeFundamentalLawSpec {
            standardized_scores_z: vec![T::one(), T::from_f64(-0.5), T::from_f64(0.7)],
            residual_volatilities: vec![T::from_f64(0.2), T::from_f64(0.25), T::from_f64(0.3)],
            benchmark_weights: vec![T::from_f64(0.3), T::from_f64(0.4), T::from_f64(0.3)],
            information_coefficient_ic: ic,
            rebalances_per_year: 12,
            signal_autocorrelation_rho: T::zero(),
            target_tracking_error: T::from_f64(0.05),
            max_active_weight_bound: T::from_f64(100.0),
        }
    }
    #[test]
    fn flam_invariants_and_dual_fd() {
        let r = evaluate_grinold_kahn_clarke_fundamental_law(&spec(0.1));
        assert!((r.autocorrelation_adjusted_breadth_br_eff - r.nominal_breadth_br).abs() < 1e-12);
        assert!((r.transfer_coefficient_tc - 1.0).abs() < 1e-12);
        assert!(
            (r.constrained_information_ratio
                - r.transfer_coefficient_tc * r.unconstrained_information_ratio)
                .abs()
                < 1e-12
        );
        let h = 1e-5;
        let fd = (evaluate_grinold_kahn_clarke_fundamental_law(&spec(0.1 + h))
            .expected_active_return
            - evaluate_grinold_kahn_clarke_fundamental_law(&spec(0.1 - h)).expected_active_return)
            / (2.0 * h);
        let ad = evaluate_grinold_kahn_clarke_fundamental_law(&spec(Dual::<1>::variable(0.1, 0)))
            .expected_active_return;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod phase27_tests {
    use super::*;
    use ql_ad::Dual;
    fn make<T: Real>(x: T) -> BussetiBoydRiskConstrainedKellySpec<T> {
        BussetiBoydRiskConstrainedKellySpec {
            excess_return_means: vec![x, T::from_f64(0.03)],
            covariance: vec![
                vec![T::from_f64(0.04), T::from_f64(0.01)],
                vec![T::from_f64(0.01), T::from_f64(0.03)],
            ],
            risk_free_rate: T::from_f64(0.02),
            drawdown_floor_alpha: T::from_f64(0.8),
            max_drawdown_prob_beta: T::from_f64(0.5),
            specialist_immediate_capital: T::from_f64(100.0),
            slow_moving_external_capital: T::from_f64(100.0),
            capital_arrival_rate: T::from_f64(0.1),
            dislocation_horizon_years: T::one(),
            forced_liquidation_notional: T::from_f64(10.0),
        }
    }
    #[test]
    fn dual_fd_parity() {
        let h = 1e-5;
        let fd = (solve_busseti_boyd_risk_constrained_kelly(&make(0.05 + h))
            .constrained_portfolio_volatility
            - solve_busseti_boyd_risk_constrained_kelly(&make(0.05 - h))
                .constrained_portfolio_volatility)
            / (2.0 * h);
        let d = solve_busseti_boyd_risk_constrained_kelly(&make(Dual::<1>::variable(0.05, 0)))
            .constrained_portfolio_volatility;
        assert!(d.d[0].abs() > 1e-6 && (d.d[0] - fd).abs() < 1e-5);
    }

    #[test]
    fn phase27_invariants() {
        let r = solve_busseti_boyd_risk_constrained_kelly(&make(0.05));
        assert!(r.kelly_fraction_kappa > 0.0 && r.kelly_fraction_kappa <= 1.0);
        assert!(r.horizon_fire_sale_discount_bps < r.initial_fire_sale_dislocation_bps);
    }
}

#[cfg(test)]
mod phase29_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn etf_and_ic_invariants_and_dual() {
        fn mk<T: Real>(x: T) -> MadhavanAbnerHillEtfApArbitrageSpec<T> {
            MadhavanAbnerHillEtfApArbitrageSpec {
                etf_market_price: x,
                reported_inav: T::from_f64(100.0),
                futures_beta_return: T::from_f64(0.01),
                stale_underlying_weight: T::from_f64(0.4),
                ap_creation_fee_bps: T::from_f64(5.0),
                basket_bid_ask_spread_bps: T::from_f64(8.0),
                creation_unit_volatility_bps: T::from_f64(10.0),
                ap_risk_aversion: T::from_f64(0.5),
                embedded_unrealized_gain_fraction: T::from_f64(0.3),
                annual_redemption_turnover: T::from_f64(0.2),
                capital_gains_tax_rate: T::from_f64(0.25),
            }
        }
        let r = evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha(&mk(101.0));
        assert!(r.annual_in_kind_heartbeat_tax_alpha_bps > 0.0);
        assert!(r.true_arbitrage_basis_bps.abs() < r.illusionary_inav_premium_bps.abs());
        let h = 1e-5;
        let fd = (evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha(&mk(101.0 + h))
            .true_arbitrage_basis_bps
            - evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha(&mk(101.0 - h))
                .true_arbitrage_basis_bps)
            / (2.0 * h);
        let ad = evaluate_madhavan_abner_hill_etf_ap_arbitrage_and_tax_alpha(&mk(
            Dual::<1>::variable(101.0, 0),
        ))
        .true_arbitrage_basis_bps;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
        let s = GrinoldQianChristophersonAlphaDecaySpec {
            fast_signal_ic0: 0.1,
            fast_signal_half_life_months: 3.0,
            slow_signal_ic0: 0.08,
            slow_signal_half_life_months: 12.0,
            signal_correlation: 0.2,
            rebalance_horizon_months: 1.0,
            transaction_cost_penalty: 0.01,
            factor_active_exposures: vec![0.1],
            factor_benchmark_returns: vec![0.02],
            normal_factor_exposures: vec![0.05],
        };
        let q = evaluate_grinold_qian_christopherson_ic_decay_and_timing(&s);
        assert!(
            (q.optimal_fast_signal_weight + q.optimal_slow_signal_weight - 1.0).abs() < 1e-12
                && (q.christopherson_policy_return + q.christopherson_factor_timing_return
                    - q.total_factor_contribution)
                    .abs()
                    < 1e-12
        );
        let sd = GrinoldQianChristophersonAlphaDecaySpec {
            fast_signal_ic0: Dual::<1>::variable(0.1, 0),
            fast_signal_half_life_months: Dual::from_value(3.0),
            slow_signal_ic0: Dual::from_value(0.08),
            slow_signal_half_life_months: Dual::from_value(12.0),
            signal_correlation: Dual::from_value(0.2),
            rebalance_horizon_months: Dual::from_value(1.0),
            transaction_cost_penalty: Dual::from_value(0.01),
            factor_active_exposures: vec![Dual::from_value(0.1)],
            factor_benchmark_returns: vec![Dual::from_value(0.02)],
            normal_factor_exposures: vec![Dual::from_value(0.05)],
        };
        let h = 1e-5;
        let fd_ic = (evaluate_grinold_qian_christopherson_ic_decay_and_timing(
            &GrinoldQianChristophersonAlphaDecaySpec {
                fast_signal_ic0: 0.1 + h,
                ..s.clone()
            },
        )
        .composite_blended_ic
            - evaluate_grinold_qian_christopherson_ic_decay_and_timing(
                &GrinoldQianChristophersonAlphaDecaySpec {
                    fast_signal_ic0: 0.1 - h,
                    ..s.clone()
                },
            )
            .composite_blended_ic)
            / (2.0 * h);
        let ad_ic =
            evaluate_grinold_qian_christopherson_ic_decay_and_timing(&sd).composite_blended_ic;
        assert!(ad_ic.d[0].abs() > 1e-6 && (ad_ic.d[0] - fd_ic).abs() < 1e-5);
    }
}
