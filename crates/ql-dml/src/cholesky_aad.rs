use ql_ad::Real;
use ql_mc::PhiloxRng;

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
    let f = s.spot * (-s.risk_free_rate * s.expiry).exp();
    let d = ((f / s.exotic_strike).ln() + T::from_f64(0.5) * sig * sig * s.expiry)
        / (sig * s.expiry.sqrt());
    let price = (-s.risk_free_rate * s.expiry).exp()
        * (f * d.norm_cdf() - s.exotic_strike * (d - sig * s.expiry.sqrt()).norm_cdf())
        * s.exotic_strike.pow(s.exotic_power_p);
    let vega = (-s.risk_free_rate * s.expiry).exp()
        * f
        * d.norm_pdf()
        * s.expiry.sqrt()
        * s.exotic_strike.pow(s.exotic_power_p);
    let g = [vega, vega * x, vega * x * x];
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
            z -= j[i][k] * lam[k];
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

#[cfg(test)]
mod phase21_tests {
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
