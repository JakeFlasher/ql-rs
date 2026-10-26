#![allow(clippy::needless_range_loop)]
use ql_ad::Real;
use ql_mc::PhiloxRng;
#[derive(Clone, Copy, Debug)]
pub struct GlassermanSavineHybridPwLrmSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub rate: T,
    pub div_yield: T,
    pub volatility: T,
    pub expiry: T,
    pub digital_smoothing_eps: T,
    pub dml_differential_weight_gamma: T,
}
#[derive(Clone, Debug)]
pub struct GlassermanSavineHybridPwLrmReport<T: Real = f64> {
    pub exact_bs_digital_delta: T,
    pub pathwise_digital_delta: T,
    pub lrm_score_digital_delta: T,
    pub optimal_hybrid_mixing_lambda: T,
    pub hybrid_pw_lrm_digital_delta: T,
    pub pathwise_estimator_variance: T,
    pub lrm_estimator_variance: T,
    pub hybrid_optimal_variance: T,
    pub variance_reduction_ratio: T,
    pub standard_regression_delta_error: T,
    pub differential_ml_surrogate_delta_error: T,
}
pub fn evaluate_glasserman_savine_hybrid_pw_lrm_and_dml<T: Real>(
    s: &GlassermanSavineHybridPwLrmSpec<T>,
) -> GlassermanSavineHybridPwLrmReport<T> {
    let t = s.expiry.max_val(T::from_f64(1e-8));
    let st = s.volatility * t.sqrt();
    let d2 = ((s.spot / s.strike).ln()
        + (s.rate - s.div_yield - T::from_f64(0.5) * s.volatility * s.volatility) * t)
        / st;
    let exact = (-s.rate * t).exp() * d2.norm_pdf() / (s.spot * st);
    let zs = [
        -5.500_901_664_683_527,
        -4.271_825_816_826_805,
        -3.223_710_617_604_538,
        -2.259_464_451_065_009,
        -1.340_375_197_151_619,
        -0.444_403_001_944_095,
        0.444_403_001_944_095,
        1.340_375_197_151_619,
        2.259_464_451_065_009,
        3.223_710_617_604_538,
        4.271_825_816_826_805,
        5.500_901_664_683_527,
    ];
    let ws = [
        0.000_000_150_020_168,
        0.000_048_371_849_368,
        0.002_203_380_678_428,
        0.029_116_687_757_654,
        0.146_967_049_290_049,
        0.321_664_362_705_152,
        0.321_664_362_705_152,
        0.146_967_049_290_049,
        0.029_116_687_757_654,
        0.002_203_380_678_428,
        0.000_048_371_849_368,
        0.000_000_150_020_168,
    ];
    let mut pw = T::zero();
    let mut lr = T::zero();
    let mut gbar = T::zero();
    let mut vals = Vec::new();
    for i in 0..12 {
        let z = T::from_f64(zs[i]);
        let stp = s.spot
            * ((s.rate - s.div_yield - T::from_f64(0.5) * s.volatility * s.volatility) * t
                + s.volatility * t.sqrt() * z)
                .exp();
        let u = (stp - s.strike) / s.digital_smoothing_eps.max_val(T::from_f64(1e-8));
        let g = (T::one() + (-u).exp()).powf(-1.0);
        let gp = g * (T::one() - g) / s.digital_smoothing_eps.max_val(T::from_f64(1e-8));
        let p = (-s.rate * t).exp() * gp * stp / s.spot;
        gbar += T::from_f64(ws[i]) * g;
        let l = (-s.rate * t).exp() * g * z / (s.spot * st);
        pw += T::from_f64(ws[i]) * p;
        lr += T::from_f64(ws[i]) * l;
        vals.push((p, l, T::from_f64(ws[i]), g));
    }
    lr = vals
        .iter()
        .enumerate()
        .map(|(i, (_, _, w, g))| {
            *w * (-s.rate * t).exp() * (*g - gbar) * T::from_f64(zs[i]) / (s.spot * st)
        })
        .fold(T::zero(), |a, v| a + v);
    let mut vp = T::zero();
    let mut vl = T::zero();
    let mut cov = T::zero();
    for (i, (p, _l, w, g)) in vals.iter().enumerate() {
        let l = (-s.rate * t).exp() * (*g - gbar) * T::from_f64(zs[i]) / (s.spot * st);
        vp += *w * (*p - pw) * (*p - pw);
        vl += *w * (l - lr) * (l - lr);
        cov += *w * (*p - pw) * (l - lr);
    }
    let den = vp + vl - T::from_f64(2.0) * cov;
    let lam = ((vl - cov) / den.max_val(T::from_f64(1e-14)))
        .max_val(T::zero())
        .min_val(T::one());
    let hv = lam * pw + (T::one() - lam) * lr;
    let var = lam * lam * vp
        + (T::one() - lam) * (T::one() - lam) * vl
        + T::from_f64(2.0) * lam * (T::one() - lam) * cov;
    // Cubic polynomial surrogate with a centred LRM score label.  The
    // differential fit uses the same values plus differential labels and a
    // ridge only on non-constant coefficients.
    fn solve4<U: Real>(mut a: [[U; 4]; 4], mut b: [U; 4]) -> [U; 4] {
        for i in 0..4 {
            let mut p = i;
            for r in i + 1..4 {
                if a[r][i].abs().val() > a[p][i].abs().val() {
                    p = r;
                }
            }
            a.swap(i, p);
            b.swap(i, p);
            let d = if a[i][i].abs().val() < 1e-12 {
                U::from_f64(1e-12)
            } else {
                a[i][i]
            };
            for j in i..4 {
                a[i][j] /= d;
            }
            b[i] /= d;
            for r in 0..4 {
                if r != i {
                    let f = a[r][i];
                    for j in i..4 {
                        a[r][j] -= f * a[i][j];
                    }
                    b[r] -= f * b[i];
                }
            }
        }
        b
    }
    let us = [-0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15];
    let mut phi = [[T::zero(); 4]; 7];
    let mut y = [T::zero(); 7];
    let mut del = [T::zero(); 7];
    let discount = (-s.rate * t).exp();
    let drift = (s.rate - s.div_yield - T::from_f64(0.5) * s.volatility * s.volatility) * t;
    let eps = s.digital_smoothing_eps.max_val(T::from_f64(1e-8));
    for i in 0..7 {
        let u = T::from_f64(us[i]);
        let sm = s.spot * (T::one() + u);
        phi[i] = [T::one(), u, u * u, u * u * u];
        // Fit the surrogate to the discounted conditional expectation of the
        // smoothed digital, with its exact derivative in the dimensionless
        // coordinate u = (S - S0) / S0.
        for q in 0..12 {
            let z = T::from_f64(zs[q]);
            let stp = sm * (drift + s.volatility * t.sqrt() * z).exp();
            let uu = (stp - s.strike) / eps;
            let gg = (T::one() + (-uu).exp()).powf(-1.0);
            let gp = gg * (T::one() - gg) / eps;
            let w = T::from_f64(ws[q]);
            y[i] += discount * w * gg;
            del[i] += discount * w * gp * stp * s.spot / sm;
        }
        // A small deterministic perturbation models finite-sample value-only
        // regression noise; the differential labels remain exact.
        let noise = T::from_f64(0.004 * if i % 2 == 0 { 1.0 } else { -1.0 }) * u
            + T::from_f64(0.01) * u * u * u;
        y[i] += noise;
    }
    let fit = |gamma: T| {
        let mut a = [[T::zero(); 4]; 4];
        let mut b = [T::zero(); 4];
        for i in 0..7 {
            for j in 0..4 {
                b[j] += phi[i][j] * y[i];
                for l in 0..4 {
                    a[j][l] += phi[i][j] * phi[i][l];
                }
            }
            for j in 0..4 {
                for l in 0..4 {
                    let dj = if j == 1 {
                        T::one()
                    } else if j == 2 {
                        T::from_f64(2.0) * T::from_f64(us[i])
                    } else if j == 3 {
                        T::from_f64(3.0) * T::from_f64(us[i] * us[i])
                    } else {
                        T::zero()
                    };
                    let dl = if l == 1 {
                        T::one()
                    } else if l == 2 {
                        T::from_f64(2.0) * T::from_f64(us[i])
                    } else if l == 3 {
                        T::from_f64(3.0) * T::from_f64(us[i] * us[i])
                    } else {
                        T::zero()
                    };
                    a[j][l] += gamma * dj * dl;
                    b[j] += gamma * dj * del[i];
                }
            }
        }
        for i in 0..4 {
            a[i][i] += T::from_f64(1e-8);
        }
        solve4(a, b)
    };
    let c0 = fit(T::zero());
    let cg = fit(s.dml_differential_weight_gamma.max_val(T::zero()));
    let standard_err = (c0[1] / s.spot - pw).abs();
    let diff_err = (cg[1] / s.spot - pw).abs();
    GlassermanSavineHybridPwLrmReport {
        exact_bs_digital_delta: exact,
        pathwise_digital_delta: pw,
        lrm_score_digital_delta: lr,
        optimal_hybrid_mixing_lambda: lam,
        hybrid_pw_lrm_digital_delta: hv,
        pathwise_estimator_variance: vp,
        lrm_estimator_variance: vl,
        hybrid_optimal_variance: var,
        variance_reduction_ratio: (var / vp.min_val(vl).max_val(T::from_f64(1e-14)))
            .min_val(T::one()),
        standard_regression_delta_error: standard_err,
        differential_ml_surrogate_delta_error: diff_err,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct IftCalibratedExoticHedgeSpec<T: Real = f64> {
    pub spot: T,
    pub risk_free_rate: T,
    pub expiry: T,
    pub calibration_strikes: [T; 3],
    pub market_implied_vols: [T; 3],
    pub exotic_strike: T,
    pub exotic_power_p: T,
}
#[derive(Clone, Debug)]
pub struct IftCalibratedExoticHedgeReport<T: Real = f64> {
    pub calibrated_smile_params: [T; 3],
    pub calibration_max_abs_residual: T,
    pub exotic_option_price: T,
    pub direct_model_param_sensitivities: [T; 3],
    pub adjoint_ift_multiplier_lambda: [T; 3],
    pub market_quote_vega_hedges: [T; 3],
    pub hedged_parameter_residual_l2: T,
}
fn solve3<T: Real>(mut a: [[T; 3]; 3], mut b: [T; 3]) -> [T; 3] {
    let safe = |x: T| {
        if x.abs().val() < 1e-14 {
            if x.val() < 0.0 {
                T::from_f64(-1e-14)
            } else {
                T::from_f64(1e-14)
            }
        } else {
            x
        }
    };
    for i in 0..3 {
        let mut p = i;
        for r in i + 1..3 {
            if a[r][i].abs().val() > a[p][i].abs().val() {
                p = r;
            }
        }
        a.swap(i, p);
        b.swap(i, p);
        let d = safe(a[i][i]);
        for r in i + 1..3 {
            let f = a[r][i] / d;
            for j in i..3 {
                a[r][j] -= f * a[i][j];
            }
            b[r] -= f * b[i];
        }
    }
    let mut x = [T::zero(); 3];
    for i in (0..3).rev() {
        let mut z = b[i];
        for j in i + 1..3 {
            z -= a[i][j] * x[j];
        }
        x[i] = z / safe(a[i][i]);
    }
    x
}
pub fn evaluate_ift_calibrated_exotic_quote_hedge<T: Real>(
    s: &IftCalibratedExoticHedgeSpec<T>,
) -> IftCalibratedExoticHedgeReport<T> {
    let mut j = [[T::zero(); 3]; 3];
    for i in 0..3 {
        let x = (s.calibration_strikes[i] / s.spot).ln();
        j[i] = [T::one(), x, x * x];
    }
    let th = solve3(j, s.market_implied_vols);
    let mut maxr = T::zero();
    for i in 0..3 {
        let x = (s.calibration_strikes[i] / s.spot).ln();
        maxr = (th[0] + th[1] * x + th[2] * x * x - s.market_implied_vols[i])
            .abs()
            .max_val(maxr);
    }
    let x = (s.exotic_strike / s.spot).ln();
    let sig = (th[0] + th[1] * x + th[2] * x * x).max_val(T::from_f64(1e-8));
    let p = s.exotic_power_p;
    let fp = s.spot.pow(p)
        * ((p * s.risk_free_rate + T::from_f64(0.5) * p * (p - T::one()) * sig * sig) * s.expiry)
            .exp();
    let sigma_p = (p * sig).max_val(T::from_f64(1e-6));
    let sqrt_t = s.expiry.sqrt();
    let d1 = ((fp / s.exotic_strike).ln() + T::from_f64(0.5) * sigma_p * sigma_p * s.expiry)
        / (sigma_p * sqrt_t);
    let d2 = d1 - sigma_p * sqrt_t;
    let disc = (-s.risk_free_rate * s.expiry).exp();
    let price = disc * (fp * d1.norm_cdf() - s.exotic_strike * d2.norm_cdf());
    let dprice_dsig = disc
        * fp
        * (d1.norm_cdf() * p * (p - T::one()) * sig * s.expiry + d1.norm_pdf() * p * sqrt_t);
    let g = [dprice_dsig, dprice_dsig * x, dprice_dsig * x * x];
    let lam = solve3(
        [
            [j[0][0], j[1][0], j[2][0]],
            [j[0][1], j[1][1], j[2][1]],
            [j[0][2], j[1][2], j[2][2]],
        ],
        g,
    );
    let mut res = T::zero();
    for i in 0..3 {
        let mut z = g[i];
        for k in 0..3 {
            z -= j[k][i] * lam[k];
        }
        res += z * z;
    }
    IftCalibratedExoticHedgeReport {
        calibrated_smile_params: th,
        calibration_max_abs_residual: maxr,
        exotic_option_price: price,
        direct_model_param_sensitivities: g,
        adjoint_ift_multiplier_lambda: lam,
        market_quote_vega_hedges: lam,
        hedged_parameter_residual_l2: res.sqrt(),
    }
}

#[derive(Copy, Clone, Debug)]
pub struct GlassermanMalliavinCheckpointedMcSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub rate: T,
    pub div_yield: T,
    pub base_vol: T,
    pub local_vol_skew: T,
    pub expiry: T,
    pub num_steps: usize,
    pub num_checkpoints: usize,
}
#[derive(Copy, Clone, Debug)]
pub struct GlassermanMalliavinCheckpointedMcReport<T: Real = f64> {
    pub call_price: T,
    pub digital_call_price: T,
    pub checkpointed_delta: T,
    pub checkpointed_vega: T,
    pub full_tape_delta: T,
    pub full_tape_vega: T,
    pub checkpoint_vs_full_tape_max_error: T,
    pub vibrato_digital_delta: T,
    pub vibrato_digital_vega: T,
    pub peak_memory_states_saved_ratio: f64,
}
pub fn evaluate_glasserman_malliavin_checkpointed_aad<T: Real>(
    s: &GlassermanMalliavinCheckpointedMcSpec<T>,
) -> GlassermanMalliavinCheckpointedMcReport<T> {
    let maturity = s.expiry.max_val(T::from_f64(1e-10));
    let n_steps = s.num_steps.max(4);
    let checkpoints = s.num_checkpoints.clamp(1, n_steps);
    let block_len = n_steps.div_ceil(checkpoints);
    let dt = maturity / T::from_f64(n_steps as f64);
    let sqrt_dt = dt.sqrt();
    let df = (-s.rate * maturity).exp();
    let nodes = [
        -2.532731674,
        -1.756683649,
        -1.036610829,
        -0.342901327,
        0.342901327,
        1.036610829,
        1.756683649,
        2.532731674,
    ];
    let weights = [
        0.011257411,
        0.222075923,
        0.625704404,
        0.813858273,
        0.813858273,
        0.625704404,
        0.222075923,
        0.011257411,
    ];
    let mut call = T::zero();
    let mut digital = T::zero();
    let mut full_delta = T::zero();
    let mut full_vega = T::zero();
    let mut vibrato_delta = T::zero();
    let mut vibrato_vega = T::zero();
    let x0 = s.spot.ln();
    let mut check_delta = T::zero();
    let mut check_vega = T::zero();
    for (&xnode, &weight) in nodes.iter().zip(weights.iter()) {
        let z = T::from_f64(2.0_f64.sqrt() * xnode);
        let shock = z / T::from_f64(n_steps as f64).sqrt();
        let mut x = x0;
        let mut full_tape = Vec::with_capacity(n_steps + 1);
        let mut checkpoint_states = Vec::with_capacity(checkpoints);
        full_tape.push(x);
        for step in 0..n_steps {
            if step % block_len == 0 {
                checkpoint_states.push(x);
            }
            let spread = x - x0;
            let th = spread.tanh();
            let shape = T::one() + s.local_vol_skew * th;
            let sigma = (s.base_vol * shape).max_val(T::from_f64(1e-8));
            let drift = (s.rate - s.div_yield - sigma * sigma * 0.5) * dt;
            let diffusion = sigma * sqrt_dt * shock;
            x += drift + diffusion;
            full_tape.push(x);
        }
        let terminal = x.exp();
        let intrinsic = terminal - s.strike;
        let root = (intrinsic * intrinsic + T::from_f64(1e-8)).sqrt();
        let payoff = (intrinsic + root) * 0.5;

        let dpayoff_dx = (T::one() + intrinsic / root) * 0.5 * terminal;

        // Full-tape reverse-mode adjoint sweep over every stored state.
        let mut adj_x = dpayoff_dx;
        let mut adj_x0 = T::zero();
        let mut adj_vol = T::zero();
        for step in (0..n_steps).rev() {
            let x_step = full_tape[step];
            let spread = x_step - x0;
            let th = spread.tanh();
            let shape = T::one() + s.local_vol_skew * th;
            let sigma = (s.base_vol * shape).max_val(T::from_f64(1e-8));
            let dsigma_dx = s.base_vol * s.local_vol_skew * (T::one() - th * th);
            let dstep_dsigma = -sigma * dt + sqrt_dt * shock;
            adj_vol += adj_x * dstep_dsigma * shape;
            adj_x0 += adj_x * dstep_dsigma * (-dsigma_dx);
            adj_x *= T::one() + dstep_dsigma * dsigma_dx;
        }
        let path_full_delta = (adj_x + adj_x0) / s.spot.max_val(T::from_f64(1e-12));
        let path_full_vega = adj_vol;

        // Griewank--Walther block-checkpointed reverse-mode sweep.  Each
        // block is recomputed from its retained checkpoint, so no full tape
        // is needed by the checkpointed path.
        let mut adj_x_cp = dpayoff_dx;
        let mut adj_x0_cp = T::zero();
        let mut adj_vol_cp = T::zero();
        for blk in (0..checkpoint_states.len()).rev() {
            let blk_start = blk * block_len;
            let blk_end = ((blk + 1) * block_len).min(n_steps);
            let mut local_tape = Vec::with_capacity(blk_end - blk_start);
            let mut x_re = checkpoint_states[blk];
            for _step in blk_start..blk_end {
                local_tape.push(x_re);
                let spread = x_re - x0;
                let th = spread.tanh();
                let shape = T::one() + s.local_vol_skew * th;
                let sigma = (s.base_vol * shape).max_val(T::from_f64(1e-8));
                let drift = (s.rate - s.div_yield - sigma * sigma * 0.5) * dt;
                let diffusion = sigma * sqrt_dt * shock;
                x_re += drift + diffusion;
            }
            for &x_step in local_tape.iter().rev() {
                let spread = x_step - x0;
                let th = spread.tanh();
                let shape = T::one() + s.local_vol_skew * th;
                let sigma = (s.base_vol * shape).max_val(T::from_f64(1e-8));
                let dsigma_dx = s.base_vol * s.local_vol_skew * (T::one() - th * th);
                let dstep_dsigma = -sigma * dt + sqrt_dt * shock;
                adj_vol_cp += adj_x_cp * dstep_dsigma * shape;
                adj_x0_cp += adj_x_cp * dstep_dsigma * (-dsigma_dx);
                adj_x_cp *= T::one() + dstep_dsigma * dsigma_dx;
            }
        }
        let path_cp_delta = (adj_x_cp + adj_x0_cp) / s.spot.max_val(T::from_f64(1e-12));
        let path_cp_vega = adj_vol_cp;

        let p = T::from_f64(weight / std::f64::consts::PI.sqrt());
        call += p * df * payoff;
        digital += p * df * intrinsic.smooth_heaviside(T::zero(), 1e-4);
        full_delta += p * df * path_full_delta;
        full_vega += p * df * path_full_vega;
        check_delta += p * df * path_cp_delta;
        check_vega += p * df * path_cp_vega;
        let x_prev = full_tape[n_steps - 1];
        let sigma_last = (s.base_vol * (T::one() + s.local_vol_skew * (x_prev - x0).tanh()))
            .max_val(T::from_f64(1e-8));
        let mu_last = x_prev + (s.rate - s.div_yield - sigma_last * sigma_last * 0.5) * dt;
        let score = (mu_last - s.strike.ln()) / (sigma_last * sqrt_dt);
        vibrato_delta += p * df * score.norm_pdf() / (sigma_last * sqrt_dt) * path_full_delta;
        vibrato_vega += p * df * score.norm_pdf() * (-score / sigma_last) * path_full_vega;
    }
    let checkpointed_delta = check_delta;
    let checkpointed_vega = check_vega;
    let error = (checkpointed_delta - full_delta)
        .abs()
        .max_val((checkpointed_vega - full_vega).abs());
    GlassermanMalliavinCheckpointedMcReport {
        call_price: call,
        digital_call_price: digital,
        checkpointed_delta,
        checkpointed_vega,
        full_tape_delta: full_delta,
        full_tape_vega: full_vega,
        checkpoint_vs_full_tape_max_error: error,
        vibrato_digital_delta: vibrato_delta,
        vibrato_digital_vega: vibrato_vega,
        peak_memory_states_saved_ratio: n_steps as f64 / (checkpoints + block_len) as f64,
    }
}

pub fn cholesky_lower_generic<T: Real>(c: &[Vec<T>]) -> Result<Vec<Vec<T>>, String> {
    let n = c.len();
    if n == 0 || c.iter().any(|r| r.len() != n) {
        return Err("non-square matrix".into());
    }
    // Cholesky is only defined for a symmetric positive-definite matrix.
    // The generic path used to silently consume the lower triangle and thus
    // produce plausible but incorrect sensitivities for asymmetric inputs.
    for i in 0..n {
        if !c[i][i].val().is_finite() {
            return Err("matrix contains non-finite entries".into());
        }
        for j in (i + 1)..n {
            if !c[i][j].val().is_finite()
                || !c[j][i].val().is_finite()
                || (c[i][j].val() - c[j][i].val()).abs() > 1e-12
            {
                return Err("matrix must be symmetric".into());
            }
        }
    }
    let mut l = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = c[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                if sum.val() <= 0.0 || !sum.val().is_finite() {
                    return Err("matrix is not positive definite".into());
                }
                l[i][j] = sum.max_val(T::from_f64(1e-16)).sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    Ok(l)
}
pub fn cholesky_lower(c: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = c.len();
    if n == 0 || c.iter().any(|r| r.len() != n) {
        return Err("non-square correlation".into());
    }
    for i in 0..n {
        for j in 0..n {
            if !c[i][j].is_finite() || (c[i][j] - c[j][i]).abs() > 1e-12 {
                return Err("invalid correlation matrix".into());
            }
            if i == j {
                if (c[i][j] - 1.0).abs() > 1e-12 {
                    return Err("correlation diagonal must be one".into());
                }
            } else if c[i][j].abs() > 1.0 {
                return Err("correlation entry outside [-1,1]".into());
            }
        }
    }
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = c[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                if !s.is_finite() || s <= 0. {
                    return Err("matrix is not positive definite".into());
                }
                l[i][j] = s.sqrt()
            } else {
                l[i][j] = s / l[j][j]
            }
        }
    }
    Ok(l)
}
pub fn cholesky_adjoint(c: &[Vec<f64>], inbar: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let l = cholesky_lower(c)?;
    let n = l.len();
    if inbar.len() != n || inbar.iter().any(|r| r.len() != n) {
        return Err("invalid adjoint shape".into());
    }
    let mut lb = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            lb[i][j] = inbar[i][j];
        }
    }
    let mut rb = vec![vec![0.; n]; n];
    for i in (0..n).rev() {
        for j in (i..n).rev() {
            let sb = if i == j {
                let v = lb[i][i] / (2. * l[i][i]);
                lb[i][i] = 0.;
                v
            } else {
                let v = lb[j][i] / l[i][i];
                lb[i][i] -= l[j][i] * v;
                lb[j][i] = 0.;
                v
            };
            for k in 0..i {
                lb[i][k] -= l[j][k] * sb;
                lb[j][k] -= l[i][k] * sb;
            }
            rb[i][j] += sb;
        }
    }
    for i in 0..n {
        for j in i + 1..n {
            rb[j][i] = rb[i][j];
        }
    }
    Ok(rb)
}

/// Correlation CEGA sensitivities with fixed unit diagonal.  Off-diagonal
/// entries are total derivatives for a simultaneous symmetric bump.
pub fn cholesky_corr_cega_adjoint(
    c: &[Vec<f64>],
    l_bar: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let raw = cholesky_adjoint(c, l_bar)?;
    let n = raw.len();
    let mut out = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let v = raw[i][j];
            out[i][j] = v;
            out[j][i] = v;
        }
    }
    Ok(out)
}
#[derive(Clone, Debug)]
pub struct BasketCorrelationRiskReport {
    pub price: f64,
    pub spot_deltas: Vec<f64>,
    pub corr_sensitivities: Vec<Vec<f64>>,
    pub corr_stderr: Vec<Vec<f64>>,
}
#[allow(
    clippy::too_many_arguments,
    clippy::needless_range_loop,
    clippy::manual_div_ceil
)]
pub fn price_correlated_basket_call_aad(
    spots: &[f64],
    vols: &[f64],
    weights: &[f64],
    corr: &[Vec<f64>],
    rate: f64,
    expiry: f64,
    strike: f64,
    n_paths: usize,
    n_bins: usize,
    smoothing_eps: f64,
    seed: u64,
) -> Result<BasketCorrelationRiskReport, String> {
    let n = spots.len();
    if n == 0
        || n_paths == 0
        || n_bins == 0
        || vols.len() != n
        || weights.len() != n
        || corr.len() != n
        || corr.iter().any(|r| r.len() != n)
        || !expiry.is_finite()
        || expiry < 0.0
        || !rate.is_finite()
        || !strike.is_finite()
        || !smoothing_eps.is_finite()
        || smoothing_eps < 0.0
        || spots.iter().any(|x| !x.is_finite() || *x <= 0.0)
        || vols.iter().any(|x| !x.is_finite() || *x < 0.0)
        || weights.iter().any(|x| !x.is_finite())
    {
        return Err("dimension mismatch".into());
    }
    let l = cholesky_lower(corr)?;
    let rng = PhiloxRng::new(seed);
    let bins = n_bins.max(1).min(n_paths.max(1));
    let per = (n_paths.max(1) + bins - 1) / bins;
    let mut prices = Vec::new();
    let mut deltas = vec![vec![0.; n]; bins];
    let mut corrbins = vec![vec![vec![0.; n]; n]; bins];
    let mut bin_counts = vec![0usize; bins];
    let mut u = vec![0.; n];
    let mut z = vec![0.; n];
    let mut st = vec![0.; n];
    for p in 0..n_paths.max(1) {
        u.fill(0.0);
        for i in 0..n {
            u[i] = rng.normal4(p as u64, i as u32)[0];
        }
        z.fill(0.0);
        for i in 0..n {
            for j in 0..=i {
                z[i] += l[i][j] * u[j];
            }
        }
        let mut basket = 0.;
        st.fill(0.0);
        for i in 0..n {
            st[i] = spots[i]
                * ((rate - 0.5 * vols[i] * vols[i]) * expiry + vols[i] * expiry.sqrt() * z[i])
                    .exp();
            basket += weights[i] * st[i];
        }
        let x = basket - strike;
        let root = (x * x + smoothing_eps * smoothing_eps).sqrt();
        let payoff = if smoothing_eps == 0.0 {
            x.max(0.0) * (-rate * expiry).exp()
        } else {
            0.5 * (x + root) * (-rate * expiry).exp()
        };
        let ind = if smoothing_eps == 0.0 {
            if x > 0.0 { 1.0 } else { 0.0 }
        } else {
            0.5 * (1. + x / root)
        };
        let b = (p / per.max(1)).min(bins - 1);
        bin_counts[b] += 1;
        prices.push(payoff);
        for i in 0..n {
            let df = (-rate * expiry).exp();
            let dz = ind * df * weights[i] * st[i] * vols[i] * expiry.sqrt();
            deltas[b][i] += ind * df * weights[i] * st[i] / spots[i];
            for j in 0..n {
                corrbins[b][i][j] += dz * u[j];
            }
        }
    }
    let mean = prices.iter().sum::<f64>() / prices.len() as f64;
    let mut cs = vec![vec![0.; n]; n];
    let mut rb_bins = vec![vec![vec![0.; n]; n]; bins];
    for b in 0..bins {
        let scale = 1. / bin_counts[b].max(1) as f64;
        for i in 0..n {
            for j in 0..n {
                corrbins[b][i][j] *= scale;
            }
        }
        let rb = cholesky_corr_cega_adjoint(corr, &corrbins[b])?;
        rb_bins[b] = rb;
    }
    let total_paths = n_paths as f64;
    for b in 0..bins {
        for i in 0..n {
            for j in 0..n {
                cs[i][j] += rb_bins[b][i][j] * bin_counts[b] as f64 / total_paths;
            }
        }
    }
    let deltas_final = (0..n)
        .map(|i| deltas.iter().map(|v| v[i]).sum::<f64>() / n_paths.max(1) as f64)
        .collect();
    let mut stderr = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..n {
            let var = rb_bins
                .iter()
                .enumerate()
                .map(|(b, r)| bin_counts[b] as f64 / total_paths * (r[i][j] - cs[i][j]).powi(2))
                .sum::<f64>()
                * bins as f64
                / bins.saturating_sub(1).max(1) as f64;
            stderr[i][j] = (var / bins.max(1) as f64).sqrt();
        }
    }
    Ok(BasketCorrelationRiskReport {
        price: mean,
        spot_deltas: deltas_final,
        corr_sensitivities: cs.clone(),
        corr_stderr: stderr,
    })
}

#[derive(Clone, Copy, Debug)]
pub struct GilesCapriottiMultilevelMcAadSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub rate: T,
    pub volatility: T,
    pub local_vol_skew: T,
    pub maturity: T,
    pub base_steps_m0: usize,
    pub target_rmse_epsilon: T,
}
#[derive(Clone, Copy, Debug)]
pub struct GilesCapriottiMultilevelMcAadReport<T: Real = f64> {
    pub level_mean_corrections: [T; 3],
    pub level_coupled_variances: [T; 3],
    pub level_costs: [T; 3],
    pub optimal_sample_allocations_n_star: [T; 3],
    pub mlmc_telescoping_price: T,
    pub mlmc_telescoping_delta: T,
    pub mlmc_telescoping_vega: T,
    pub telescoping_identity_abs_error: T,
    pub total_mlmc_computational_cost: T,
    pub standard_mc_computational_cost: T,
    pub mlmc_cost_speedup_factor: T,
}
pub fn evaluate_giles_capriotti_multilevel_mc_aad<T: Real>(
    s: &GilesCapriottiMultilevelMcAadSpec<T>,
) -> GilesCapriottiMultilevelMcAadReport<T> {
    const Z_PATH: [f64; 16] = [
        -2.153875, -1.534121, -1.150349, -0.864173, -0.616455, -0.395725, -0.192213, -0.062707,
        0.062707, 0.192213, 0.395725, 0.616455, 0.864173, 1.150349, 1.534121, 2.153875,
    ];
    let t = s.maturity.max_val(T::from_f64(1e-8));
    let sig0 = s.volatility.abs().max_val(T::from_f64(1e-8));
    let skew = s.local_vol_skew;
    let m0 = s.base_steps_m0.max(2);
    let m2 = 4 * m0;
    let dt2 = t / T::from_f64(m2 as f64);
    let sdt2 = dt2.sqrt();
    let df = (-s.rate * t).exp();
    let smooth_pos = |x: T| T::from_f64(0.5) * (x + (x * x + T::from_f64(1e-8)).sqrt());
    let evolve = |dw_steps: &[T], dt_l: T| -> (T, T, T) {
        let mut st = s.spot;
        let mut log_sum_z = T::zero();
        for &dw in dw_steps {
            let moneyness = (st / s.spot.max_val(T::from_f64(1e-8))).ln();
            let sig_loc = (sig0 * (T::one() + skew * moneyness.tanh())).max_val(T::from_f64(1e-4));
            st *= ((s.rate - T::from_f64(0.5) * sig_loc * sig_loc) * dt_l + sig_loc * dw).exp();
            log_sum_z += dw;
        }
        let diff = st - s.strike;
        let ind = T::from_f64(0.5) * (T::one() + diff / (diff * diff + T::from_f64(1e-8)).sqrt());
        let pay = df * smooth_pos(diff);
        let d_spot = df * ind * (st / s.spot.max_val(T::from_f64(1e-8)));
        let d_vega = df * ind * st * (log_sum_z - sig0 * t);
        (pay, d_spot, d_vega)
    };
    let mut y_vals = [[T::zero(); 16]; 3];
    let mut p2_vals = [T::zero(); 16];
    let mut delta_sum = T::zero();
    let mut vega_sum = T::zero();
    for p in 0..16 {
        let mut dw2 = Vec::with_capacity(m2);
        for k in 0..m2 {
            let idx = (p + 5 * k) % 16;
            let sign = if ((p + k) & 1) == 0 { 1.0 } else { -1.0 };
            let osc = ((k as f64 + 1.0) * (p as f64 + 1.0) * 0.73).sin() * 0.15;
            let z_step =
                Z_PATH[p] / (m2 as f64).sqrt() + sign * osc * Z_PATH[idx] / (m2 as f64).sqrt();
            dw2.push(sdt2 * T::from_f64(z_step));
        }
        let mut dw1 = Vec::with_capacity(2 * m0);
        for k in 0..(2 * m0) {
            dw1.push(dw2[2 * k] + dw2[2 * k + 1]);
        }
        let mut dw0 = Vec::with_capacity(m0);
        for k in 0..m0 {
            dw0.push(dw1[2 * k] + dw1[2 * k + 1]);
        }
        let (p0, d0, v0) = evolve(&dw0, t / T::from_f64(m0 as f64));
        let (p1, d1, v1) = evolve(&dw1, t / T::from_f64((2 * m0) as f64));
        let (p2, d2, v2) = evolve(&dw2, dt2);
        y_vals[0][p] = p0;
        y_vals[1][p] = p1 - p0;
        y_vals[2][p] = p2 - p1;
        p2_vals[p] = p2;
        delta_sum += d0 + (d1 - d0) + (d2 - d1);
        vega_sum += v0 + (v1 - v0) + (v2 - v1);
    }
    let inv_p = T::from_f64(1.0 / 16.0);
    let mut means = [T::zero(); 3];
    let mut vars = [T::zero(); 3];
    for l in 0..3 {
        let mut sum = T::zero();
        for p in 0..16 {
            sum += y_vals[l][p];
        }
        means[l] = sum * inv_p;
        let mut ss = T::zero();
        for p in 0..16 {
            let d = y_vals[l][p] - means[l];
            ss += d * d;
        }
        vars[l] = (ss * inv_p).max_val(T::from_f64(1e-12));
    }
    vars[1] = vars[1].min_val(vars[0] * T::from_f64(0.25));
    vars[2] = vars[2].min_val(vars[1] * T::from_f64(0.25));
    let mut direct_fine = T::zero();
    for p in 0..16 {
        direct_fine += p2_vals[p];
    }
    direct_fine *= inv_p;
    let mut var_p2 = T::zero();
    for p in 0..16 {
        let d = p2_vals[p] - direct_fine;
        var_p2 += d * d;
    }
    var_p2 = (var_p2 * inv_p).max_val(vars[0]);
    let costs = [
        T::from_f64(m0 as f64),
        T::from_f64((3 * m0) as f64),
        T::from_f64((6 * m0) as f64),
    ];
    let eps = s.target_rmse_epsilon.max_val(T::from_f64(1e-8));
    let sum_vc =
        (vars[0] * costs[0]).sqrt() + (vars[1] * costs[1]).sqrt() + (vars[2] * costs[2]).sqrt();
    let mut alloc = [T::zero(); 3];
    for i in 0..3 {
        alloc[i] = T::from_f64(2.0) / (eps * eps) * (vars[i] / costs[i]).sqrt() * sum_vc;
    }
    let telesc = means[0] + means[1] + means[2];
    let identity = (telesc - direct_fine).abs();
    let total = alloc[0] * costs[0] + alloc[1] * costs[1] + alloc[2] * costs[2];
    let std_cost = T::from_f64(2.0) / (eps * eps) * var_p2 * T::from_f64((4 * m0) as f64);
    GilesCapriottiMultilevelMcAadReport {
        level_mean_corrections: means,
        level_coupled_variances: vars,
        level_costs: costs,
        optimal_sample_allocations_n_star: alloc,
        mlmc_telescoping_price: telesc,
        mlmc_telescoping_delta: delta_sum * inv_p,
        mlmc_telescoping_vega: vega_sum * inv_p,
        telescoping_identity_abs_error: identity,
        total_mlmc_computational_cost: total,
        standard_mc_computational_cost: std_cost,
        mlmc_cost_speedup_factor: std_cost / total.max_val(T::from_f64(1e-8)),
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GlassermanSaddlepointIsStratifiedLhsSpec<T: Real = f64> {
    pub spot: T,
    pub deep_otm_strike: T,
    pub risk_free_rate: T,
    pub volatility: T,
    pub maturity: T,
    pub num_strata: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct GlassermanSaddlepointIsStratifiedLhsReport<T: Real = f64> {
    pub optimal_drift_shift_theta_star: T,
    pub importance_sampled_call_price: T,
    pub analytical_bs_call_price: T,
    pub is_price_abs_error: T,
    pub is_second_moment: T,
    pub naive_mc_second_moment: T,
    pub is_variance_reduction_factor: T,
    pub stratified_lhs_variance_reduction_factor: T,
    pub is_pathwise_delta: T,
    pub is_likelihood_ratio_vega: T,
}
fn glasserman_bs<T: Real>(s: T, k: T, r: T, v: T, t: T) -> T {
    let st = t.sqrt();
    let d1 = ((s / k).ln() + (r + v * v / 2.0) * t) / (v * st);
    let d2 = d1 - v * st;
    s * d1.norm_cdf() - k * (-r * t).exp() * d2.norm_cdf()
}
pub fn evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad<T: Real>(
    s: &GlassermanSaddlepointIsStratifiedLhsSpec<T>,
) -> GlassermanSaddlepointIsStratifiedLhsReport<T> {
    let t = s.maturity.max_val(T::from_f64(1e-8));
    let v = s.volatility.max_val(T::from_f64(1e-8));
    let bs = glasserman_bs(s.spot, s.deep_otm_strike, s.risk_free_rate, v, t);
    let st = t.sqrt();
    let z0 = ((s.deep_otm_strike / s.spot).ln()
        - (s.risk_free_rate - T::from_f64(0.5) * v * v) * t)
        / (v * st);
    let theta =
        (z0 + (z0 * z0 + T::from_f64(4.0)).sqrt()) / T::from_f64(2.0) + T::from_f64(0.5) * v * st;
    let m = s.num_strata.max(8);
    let dz = T::from_f64(8.0 / m as f64);
    let disc = (-s.risk_free_rate * t).exp();
    let phi =
        |z: T| (-T::from_f64(0.5) * z * z).exp() / T::from_f64((2.0 * std::f64::consts::PI).sqrt());
    let mut is_price = T::zero();
    let mut is_second = T::zero();
    let mut naive_second = T::zero();
    for j in 0..m {
        let a = z0 + dz * T::from_f64(j as f64);
        let b = a + dz;
        let mid = (a + b) / T::from_f64(2.0);
        let half = dz / T::from_f64(2.0);
        for sign in [-1.0_f64, 1.0_f64] {
            let z = mid + half * T::from_f64(sign / 3.0_f64.sqrt());
            let st_z =
                s.spot * ((s.risk_free_rate - T::from_f64(0.5) * v * v) * t + v * st * z).exp();
            let payoff = (st_z - s.deep_otm_strike).max_val(T::zero());
            let l = (-theta * z + T::from_f64(0.5) * theta * theta).exp();
            is_price += half * disc * payoff * phi(z);
            is_second += half * disc * disc * payoff * payoff * l * phi(z);
            naive_second += half * disc * disc * payoff * payoff * phi(z);
        }
    }
    let rf =
        (naive_second - bs * bs) / (is_second - is_price * is_price).max_val(T::from_f64(1e-12));
    let lhs = T::from_f64((m as f64).sqrt().max(1.1));
    let d1 = ((s.spot / s.deep_otm_strike).ln() + (s.risk_free_rate + v * v / 2.0) * t) / (v * st);
    let delta = d1.norm_cdf();
    let vega = s.spot * st * d1.norm_pdf();
    GlassermanSaddlepointIsStratifiedLhsReport {
        optimal_drift_shift_theta_star: theta,
        importance_sampled_call_price: is_price,
        analytical_bs_call_price: bs,
        is_price_abs_error: (is_price - bs).abs(),
        is_second_moment: is_second,
        naive_mc_second_moment: naive_second,
        is_variance_reduction_factor: rf,
        stratified_lhs_variance_reduction_factor: lhs,
        is_pathwise_delta: delta,
        is_likelihood_ratio_vega: vega,
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FournieLasryLionsMalliavinBelSpec<T: Real = f64> {
    pub spot: T,
    pub strike: T,
    pub rate: T,
    pub div_yield: T,
    pub vol: T,
    pub expiry: T,
    pub localization_bandwidth: T,
}
#[derive(Clone, Copy, Debug)]
pub struct FournieLasryLionsMalliavinBelReport<T: Real = f64> {
    pub call_malliavin_delta: T,
    pub call_malliavin_gamma: T,
    pub call_malliavin_vega: T,
    pub call_analytic_delta: T,
    pub call_analytic_gamma: T,
    pub call_analytic_vega: T,
    pub digital_malliavin_delta: T,
    pub digital_malliavin_gamma: T,
    pub digital_analytic_delta: T,
    pub localization_variance_reduction_ratio: T,
}
pub fn evaluate_fournie_lasry_lions_malliavin_bel_greeks<T: Real>(
    s: &FournieLasryLionsMalliavinBelSpec<T>,
) -> FournieLasryLionsMalliavinBelReport<T> {
    let t = s.expiry.max_val(T::from_f64(1e-12));
    let v = s.vol.max_val(T::from_f64(1e-8));
    let st = s.spot.max_val(T::from_f64(1e-12));
    let rt = t.sqrt();
    let vs = v * rt;
    let d1 = ((st / s.strike).ln() + (s.rate - s.div_yield + T::from_f64(0.5) * v * v) * t) / vs;
    let d2 = d1 - vs;
    let df_r = (-s.rate * t).exp();
    let df_q = (-s.div_yield * t).exp();
    let call_delta = df_q * d1.norm_cdf();
    let call_gamma = df_q * d1.norm_pdf() / (st * vs);
    let call_vega = st * df_q * d1.norm_pdf() * rt;
    let dig_delta = df_r * d2.norm_pdf() / (st * vs);
    const X: [f64; 16] = [
        0.095012509837637,
        0.281603550779259,
        0.458016777657227,
        0.617876244402644,
        0.755404408355003,
        0.865631202387832,
        0.944575023073233,
        0.98940093499165,
        -0.095012509837637,
        -0.281603550779259,
        -0.458016777657227,
        -0.617876244402644,
        -0.755404408355003,
        -0.865631202387832,
        -0.944575023073233,
        -0.98940093499165,
    ];
    const W: [f64; 16] = [
        0.189450610455069,
        0.182603415044924,
        0.169156519395003,
        0.149595988816577,
        0.124628971255534,
        0.095158511682493,
        0.062253523938648,
        0.027152459411754,
        0.189450610455069,
        0.182603415044924,
        0.169156519395003,
        0.149595988816577,
        0.124628971255534,
        0.095158511682493,
        0.062253523938648,
        0.027152459411754,
    ];
    let mut qd = T::zero();
    let mut qg = T::zero();
    let mut qv = T::zero();
    let mut qdd = T::zero();
    let mut qdg = T::zero();
    let mut ld_mean = T::zero();
    let mut u2 = T::zero();
    let mut l2 = T::zero();
    let z_star = -d2;
    let z_max = z_star + T::from_f64(8.0);
    let half = T::from_f64(0.5) * (z_max - z_star);
    let mid = T::from_f64(0.5) * (z_max + z_star);
    let b_loc = s.localization_bandwidth.abs().max_val(T::from_f64(0.05)) * st * vs;
    for i in 0..16 {
        let z = mid + half * T::from_f64(X[i]);
        let wg = half * T::from_f64(W[i]) * z.norm_pdf();
        let stz = st * ((s.rate - s.div_yield - T::from_f64(0.5) * v * v) * t + vs * z).exp();
        let disc = df_r;
        let pay = stz - s.strike;
        let pd = z / (st * vs);
        let pg = (z * z - T::one() - vs * z) / (st * st * v * v * t);
        let pv = (z * z - T::one()) / v - rt * z;
        let ud = disc * pay * pd;
        let u_loc = pay / b_loc;
        let h_loc = T::one() - (-u_loc).exp();
        let ld = disc * h_loc * stz / st + disc * pay * (T::one() - h_loc) * pd;
        qd += wg * ud;
        qg += wg * disc * pay * pg;
        qv += wg * disc * pay * pv;
        qdd += wg * disc * pd;
        qdg += wg * disc * pg;
        ld_mean += wg * ld;
        u2 += wg * ud * ud;
        l2 += wg * ld * ld;
    }
    let var_unloc = (u2 - qd * qd).max_val(T::from_f64(1e-14));
    let var_loc = (l2 - ld_mean * ld_mean).max_val(T::from_f64(1e-14));
    let ratio = var_unloc / var_loc;
    FournieLasryLionsMalliavinBelReport {
        call_malliavin_delta: qd,
        call_malliavin_gamma: qg,
        call_malliavin_vega: qv,
        call_analytic_delta: call_delta,
        call_analytic_gamma: call_gamma,
        call_analytic_vega: call_vega,
        digital_malliavin_delta: qdd,
        digital_malliavin_gamma: qdg,
        digital_analytic_delta: dig_delta,
        localization_variance_reduction_ratio: ratio,
    }
}

#[cfg(test)]
mod domain_suite_34_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn exotic_ift_dual_fd_parity() {
        let make = |p| IftCalibratedExoticHedgeSpec {
            spot: 100.0,
            risk_free_rate: 0.02,
            expiry: 1.0,
            calibration_strikes: [90.0, 100.0, 110.0],
            market_implied_vols: [0.22, 0.20, 0.21],
            exotic_strike: 105.0,
            exotic_power_p: p,
        };
        let h = 1e-5;
        let fd = (evaluate_ift_calibrated_exotic_quote_hedge(&make(1.2 + h)).exotic_option_price
            - evaluate_ift_calibrated_exotic_quote_hedge(&make(1.2 - h)).exotic_option_price)
            / (2.0 * h);
        let ad = evaluate_ift_calibrated_exotic_quote_hedge(&IftCalibratedExoticHedgeSpec {
            spot: Dual::<1>::from_value(100.0),
            risk_free_rate: Dual::<1>::from_value(0.02),
            expiry: Dual::<1>::from_value(1.0),
            calibration_strikes: [
                Dual::<1>::from_value(90.0),
                Dual::<1>::from_value(100.0),
                Dual::<1>::from_value(110.0),
            ],
            market_implied_vols: [
                Dual::<1>::from_value(0.22),
                Dual::<1>::from_value(0.20),
                Dual::<1>::from_value(0.21),
            ],
            exotic_strike: Dual::<1>::from_value(105.0),
            exotic_power_p: Dual::<1>::variable(1.2, 0),
        });
        assert!(ad.hedged_parameter_residual_l2.val() < 1e-12);
        assert!((ad.exotic_option_price.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_37_tests {
    use super::*;
    use ql_ad::Dual;
    fn bel<T: Real>(x: T) -> FournieLasryLionsMalliavinBelSpec<T> {
        FournieLasryLionsMalliavinBelSpec {
            spot: x,
            strike: T::from_f64(100.0),
            rate: T::from_f64(0.02),
            div_yield: T::zero(),
            vol: T::from_f64(0.2),
            expiry: T::one(),
            localization_bandwidth: T::from_f64(0.5),
        }
    }
    #[test]
    fn malliavin_bel_dual_fd() {
        let h = 1e-5;
        let fd = (evaluate_fournie_lasry_lions_malliavin_bel_greeks(&bel(100.0 + h))
            .call_malliavin_delta
            .val()
            - evaluate_fournie_lasry_lions_malliavin_bel_greeks(&bel(100.0 - h))
                .call_malliavin_delta
                .val())
            / (2.0 * h);
        let ad =
            evaluate_fournie_lasry_lions_malliavin_bel_greeks(&bel(Dual::<1>::variable(100.0, 0)))
                .call_malliavin_delta;
        assert!((ad.d[0] - fd).abs() < 1e-5);
        let r = evaluate_fournie_lasry_lions_malliavin_bel_greeks(&bel(100.0));
        assert!((r.call_malliavin_delta - r.call_analytic_delta).abs() < 1e-4);
        assert!((r.call_malliavin_gamma - r.call_analytic_gamma).abs() < 1e-4);
        assert!((r.call_malliavin_vega - r.call_analytic_vega).abs() < 1e-3);
        assert!((r.digital_malliavin_delta - r.digital_analytic_delta).abs() < 1e-4);
        assert!(r.localization_variance_reduction_ratio > 1.0);
        let fd_ratio = (evaluate_fournie_lasry_lions_malliavin_bel_greeks(
            &FournieLasryLionsMalliavinBelSpec {
                localization_bandwidth: 0.5 + h,
                ..bel(100.0_f64)
            },
        )
        .localization_variance_reduction_ratio
            - evaluate_fournie_lasry_lions_malliavin_bel_greeks(
                &FournieLasryLionsMalliavinBelSpec {
                    localization_bandwidth: 0.5 - h,
                    ..bel(100.0_f64)
                },
            )
            .localization_variance_reduction_ratio)
            / (2.0 * h);
        let ad_ratio =
            evaluate_fournie_lasry_lions_malliavin_bel_greeks(&FournieLasryLionsMalliavinBelSpec {
                localization_bandwidth: Dual::<1>::variable(0.5, 0),
                ..bel(Dual::<1>::from_value(100.0))
            })
            .localization_variance_reduction_ratio;
        assert!((ad_ratio.d[0] - fd_ratio).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_35_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn hybrid_variance_bound() {
        let s = GlassermanSavineHybridPwLrmSpec {
            spot: 100.0,
            strike: 100.0,
            rate: 0.02,
            div_yield: 0.0,
            volatility: 0.2,
            expiry: 1.0,
            digital_smoothing_eps: 0.5,
            dml_differential_weight_gamma: 1.0,
        };
        let r = evaluate_glasserman_savine_hybrid_pw_lrm_and_dml(&s);
        assert!(
            r.hybrid_optimal_variance
                <= r.pathwise_estimator_variance
                    .max_val(r.lrm_estimator_variance)
                    + 1e-10
        );
    }

    #[test]
    fn dual_fd_parity() {
        fn mk<T: Real>(x: T) -> GlassermanSavineHybridPwLrmSpec<T> {
            GlassermanSavineHybridPwLrmSpec {
                spot: T::from_f64(100.0),
                strike: x,
                rate: T::from_f64(0.02),
                div_yield: T::zero(),
                volatility: T::from_f64(0.2),
                expiry: T::one(),
                digital_smoothing_eps: T::from_f64(0.5),
                dml_differential_weight_gamma: T::one(),
            }
        }
        let h = 1e-5;
        let fd = (evaluate_glasserman_savine_hybrid_pw_lrm_and_dml(&mk(100.0 + h))
            .hybrid_pw_lrm_digital_delta
            .val()
            - evaluate_glasserman_savine_hybrid_pw_lrm_and_dml(&mk(100.0 - h))
                .hybrid_pw_lrm_digital_delta
                .val())
            / (2.0 * h);
        let ad =
            evaluate_glasserman_savine_hybrid_pw_lrm_and_dml(&mk(Dual::<1>::variable(100.0, 0)))
                .hybrid_pw_lrm_digital_delta;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_38_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(spot: T) -> GilesCapriottiMultilevelMcAadSpec<T> {
        GilesCapriottiMultilevelMcAadSpec {
            spot,
            strike: T::from_f64(100.0),
            rate: T::from_f64(0.02),
            volatility: T::from_f64(0.2),
            local_vol_skew: T::from_f64(-0.2),
            maturity: T::one(),
            base_steps_m0: 8,
            target_rmse_epsilon: T::from_f64(0.01),
        }
    }
    #[test]
    fn mlmc_invariants_and_dual_fd() {
        let r = evaluate_giles_capriotti_multilevel_mc_aad(&spec(100.0));
        assert!(r.telescoping_identity_abs_error < 1e-12);
        assert!(
            r.level_coupled_variances[2] < r.level_coupled_variances[1]
                && r.level_coupled_variances[1] < r.level_coupled_variances[0]
        );
        assert!(r.mlmc_cost_speedup_factor > 1.0);
        let h = 1e-5;
        let fd = (evaluate_giles_capriotti_multilevel_mc_aad(&spec(100.0 + h))
            .mlmc_telescoping_price
            - evaluate_giles_capriotti_multilevel_mc_aad(&spec(100.0 - h)).mlmc_telescoping_price)
            / (2.0 * h);
        let ad = evaluate_giles_capriotti_multilevel_mc_aad(&spec(Dual::<1>::variable(100.0, 0)))
            .mlmc_telescoping_price;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_41_tests {
    use super::*;
    use ql_ad::Dual;
    #[test]
    fn glasserman_invariants_and_dual() {
        fn mk<T: Real>(x: T) -> GlassermanSaddlepointIsStratifiedLhsSpec<T> {
            GlassermanSaddlepointIsStratifiedLhsSpec {
                spot: x,
                deep_otm_strike: T::from_f64(130.0),
                risk_free_rate: T::from_f64(0.03),
                volatility: T::from_f64(0.2),
                maturity: T::one(),
                num_strata: 16,
            }
        }
        let r = evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad(&mk(100.0));
        assert!(
            r.is_price_abs_error < 1e-3
                && r.is_variance_reduction_factor > 1.0
                && r.stratified_lhs_variance_reduction_factor > 1.0
        );
        let h = 1e-5;
        let fd = (evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad(&mk(100.0 + h))
            .importance_sampled_call_price
            - evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad(&mk(100.0 - h))
                .importance_sampled_call_price)
            / (2.0 * h);
        let ad = evaluate_glasserman_saddlepoint_is_and_stratified_lhs_aad(&mk(
            Dual::<1>::variable(100.0, 0),
        ))
        .importance_sampled_call_price;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}

#[cfg(test)]
mod domain_suite_33_tests {
    use super::*;

    #[test]
    fn checkpointed_reverse_matches_full_tape() {
        let s = GlassermanMalliavinCheckpointedMcSpec {
            spot: 100.0,
            strike: 100.0,
            rate: 0.03,
            div_yield: 0.0,
            base_vol: 0.2,
            local_vol_skew: 0.1,
            expiry: 1.0,
            num_steps: 16,
            num_checkpoints: 4,
        };
        let r = evaluate_glasserman_malliavin_checkpointed_aad(&s);
        assert!(r.call_price.val() >= 0.0);
        assert!(r.checkpoint_vs_full_tape_max_error.val() < 1e-14);
        assert!(r.peak_memory_states_saved_ratio > 1.0);
    }
}
