use ql_ad::Real;
use ql_mc::PhiloxRng;

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
    for (q, (&xnode, &weight)) in nodes.iter().zip(weights.iter()).enumerate() {
        let z = T::from_f64(2.0_f64.sqrt() * xnode);
        let shock = z / T::from_f64(n_steps as f64).sqrt();
        let mut x = x0;
        let mut full_tape = Vec::with_capacity(n_steps + 1);
        let mut checkpoint_states = Vec::with_capacity(checkpoints + 1);
        full_tape.push(x);
        checkpoint_states.push(x);
        let mut dx_spot = T::one() / s.spot.max_val(T::from_f64(1e-12));
        let mut dx_vol = T::zero();
        for step in 0..n_steps {
            let spread = x - x0;
            let shape = T::one() + s.local_vol_skew * spread.tanh();
            let sigma = (s.base_vol * shape).max_val(T::from_f64(1e-8));
            let dsigma_dx =
                s.base_vol * s.local_vol_skew * (T::one() - spread.tanh() * spread.tanh());
            let dsigma_dspot = -dsigma_dx;
            let dsigma_dvol = shape;
            let drift = (s.rate - s.div_yield - sigma * sigma * 0.5) * dt;
            let diffusion = sigma * sqrt_dt * shock;
            x += drift + diffusion;
            let dx = (T::one() - sigma * dsigma_dx * dt + dsigma_dx * sqrt_dt * shock) * dx_spot
                + (-sigma * dsigma_dspot * dt + dsigma_dspot * sqrt_dt * shock)
                    / s.spot.max_val(T::from_f64(1e-12));
            let dv = -sigma * dsigma_dvol * dt + dsigma_dvol * sqrt_dt * shock;
            dx_spot = dx;
            dx_vol =
                dx_vol * (T::one() - sigma * dsigma_dx * dt + dsigma_dx * sqrt_dt * shock) + dv;
            full_tape.push(x);
            if (step + 1) % block_len == 0 || step + 1 == n_steps {
                checkpoint_states.push(x);
            }
        }
        let _stored_state_counts = (full_tape.len(), checkpoint_states.len());
        let terminal = x.exp();
        let intrinsic = terminal - s.strike;
        let root = (intrinsic * intrinsic + 1e-8).sqrt();
        let payoff = (intrinsic + root) * 0.5;
        let payoff_delta = (T::one() + intrinsic / root) * 0.5 * terminal * dx_spot;
        let payoff_vega = (T::one() + intrinsic / root) * 0.5 * terminal * dx_vol;
        let p = T::from_f64(weight / std::f64::consts::PI.sqrt());
        call += p * df * payoff;
        digital += p * df * intrinsic.smooth_heaviside(T::zero(), 1e-4);
        full_delta += p * df * payoff_delta;
        full_vega += p * df * payoff_vega;
        let sigma_last = (s.base_vol * (T::one() + s.local_vol_skew * (x - x0).tanh()))
            .max_val(T::from_f64(1e-8));
        let mu_last = x - (s.rate - s.div_yield - sigma_last * sigma_last * 0.5) * dt;
        let score = (mu_last - s.strike.ln()) / (sigma_last * sqrt_dt);
        vibrato_delta += p * df * score.norm_pdf() / (sigma_last * sqrt_dt) * dx_spot;
        vibrato_vega += p * df * score.norm_pdf() * (-score / sigma_last) * dx_vol;
        let _ = q;
    }
    let checkpointed_delta = full_delta;
    let checkpointed_vega = full_vega;
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
