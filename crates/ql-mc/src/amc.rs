use crate::PhiloxRng;
use nalgebra::{DMatrix, DVector};
use ql_ad::Real;
#[derive(Clone, Debug)]
pub struct ExerciseBoundaryPolicy {
    pub coeffs: Vec<Vec<f64>>,
    pub spot_scale: f64,
}
#[derive(Copy, Clone, Debug)]
pub struct AmcEngine {
    pub n_calibration_paths: usize,
    pub n_valuation_paths: usize,
    pub n_steps: usize,
    pub basis_degree: usize,
    pub seed: u64,
}
fn basis(x: f64, d: usize) -> Vec<f64> {
    let mut b = vec![0.; d + 1];
    b[0] = 1.;
    if d > 0 {
        b[1] = x;
    }
    for k in 2..=d {
        b[k] = 2. * x * b[k - 1] - b[k - 2];
    }
    b
}
fn intrinsic(s: f64, k: f64, call: bool) -> f64 {
    if call {
        (s - k).max(0.)
    } else {
        (k - s).max(0.)
    }
}
impl AmcEngine {
    #[allow(clippy::needless_range_loop)]
    pub fn calibrate_policy(
        &self,
        strike: f64,
        spot: f64,
        vol: f64,
        rate: f64,
        expiry: f64,
        is_call: bool,
    ) -> ExerciseBoundaryPolicy {
        self.calibrate_policy_with_div(strike, spot, vol, rate, 0.0, expiry, is_call)
    }

    /// Calibrate a Longstaff--Schwartz policy under a continuous dividend
    /// yield.  The policy basis is unchanged; only the risk-neutral drift is
    /// shifted from `rate` to `rate-div`.
    #[allow(clippy::needless_range_loop)]
    #[allow(clippy::too_many_arguments)]
    pub fn calibrate_policy_with_div(
        &self,
        strike: f64,
        spot: f64,
        vol: f64,
        rate: f64,
        div: f64,
        expiry: f64,
        is_call: bool,
    ) -> ExerciseBoundaryPolicy {
        let n = self.n_steps.max(1);
        let paths = self.n_calibration_paths.max(self.basis_degree + 2);
        let dt = expiry / n as f64;
        let mut s = vec![vec![0.; n + 1]; paths];
        let rng = PhiloxRng::new(self.seed);
        for p in 0..paths {
            s[p][0] = spot;
            for j in 0..n {
                let z = rng.normal4(p as u64, j as u32)[0];
                s[p][j + 1] =
                    s[p][j] * ((rate - div - 0.5 * vol * vol) * dt + vol * dt.sqrt() * z).exp();
            }
        }
        let mut cf = (0..paths)
            .map(|p| intrinsic(s[p][n], strike, is_call))
            .collect::<Vec<_>>();
        let mut coeffs = vec![Vec::new(); n];
        for j in (1..n).rev() {
            for x in &mut cf {
                *x *= (-rate * dt).exp();
            }
            let itm = (0..paths)
                .filter(|&p| intrinsic(s[p][j], strike, is_call) > 0.)
                .collect::<Vec<_>>();
            if itm.len() > self.basis_degree + 1 {
                let mut a = DMatrix::zeros(itm.len(), self.basis_degree + 1);
                let mut y = DVector::zeros(itm.len());
                for (r, &p) in itm.iter().enumerate() {
                    let b = basis(s[p][j] / strike - 1., self.basis_degree);
                    for c in 0..b.len() {
                        a[(r, c)] = b[c];
                    }
                    y[r] = cf[p];
                }
                if let Ok(beta) = a.svd(true, true).solve(&y, 1e-10) {
                    coeffs[j] = beta.iter().copied().collect();
                    for &p in &itm {
                        let b = basis(s[p][j] / strike - 1., self.basis_degree);
                        let cont = coeffs[j]
                            .iter()
                            .zip(b.iter())
                            .map(|(x, y)| x * y)
                            .sum::<f64>();
                        let g = intrinsic(s[p][j], strike, is_call);
                        if g >= cont {
                            cf[p] = g;
                        }
                    }
                }
            }
        }
        ExerciseBoundaryPolicy {
            coeffs,
            spot_scale: spot.max(strike),
        }
    }

    pub fn calibrate_put_policy_with_div(
        &self,
        strike: f64,
        spot: f64,
        vol: f64,
        rate: f64,
        div: f64,
        expiry: f64,
    ) -> ExerciseBoundaryPolicy {
        self.calibrate_policy_with_div(strike, spot, vol, rate, div, expiry, false)
    }
    #[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
    pub fn price_with_policy(
        &self,
        policy: &ExerciseBoundaryPolicy,
        strike: f64,
        spot: f64,
        vol: f64,
        rate: f64,
        expiry: f64,
        is_call: bool,
    ) -> (f64, f64) {
        self.price_with_policy_and_div(policy, strike, spot, vol, rate, 0.0, expiry, is_call)
    }

    #[allow(clippy::too_many_arguments, clippy::needless_range_loop)]
    pub fn price_with_policy_and_div(
        &self,
        policy: &ExerciseBoundaryPolicy,
        strike: f64,
        spot: f64,
        vol: f64,
        rate: f64,
        div: f64,
        expiry: f64,
        is_call: bool,
    ) -> (f64, f64) {
        let n = self.n_steps.max(1);
        let dt = expiry / n as f64;
        let rng = PhiloxRng::new(self.seed + 1_000_000);
        let mut sum = 0.;
        let mut ss = 0.;
        for p in 0..self.n_valuation_paths.max(1) {
            let mut s = spot;
            let mut payoff = None;
            for j in 1..=n {
                let z = rng.normal4(p as u64, j as u32)[0];
                s *= ((rate - div - 0.5 * vol * vol) * dt + vol * dt.sqrt() * z).exp();
                let g = intrinsic(s, strike, is_call);
                if j < n && !policy.coeffs[j].is_empty() && g > 0. {
                    let b = basis(s / strike - 1., self.basis_degree);
                    let c = policy.coeffs[j]
                        .iter()
                        .zip(b.iter())
                        .map(|(x, y)| x * y)
                        .sum::<f64>();
                    if g >= c {
                        payoff = Some(g * (-rate * j as f64 * dt).exp());
                        break;
                    }
                }
                if j == n {
                    payoff = Some(g * (-rate * expiry).exp());
                }
            }
            let x = payoff.unwrap_or(0.);
            sum += x;
            ss += x * x;
        }
        let mean = sum / self.n_valuation_paths.max(1) as f64;
        let se = ((ss / self.n_valuation_paths.max(1) as f64 - mean * mean).max(0.)
            / self.n_valuation_paths.max(1) as f64)
            .sqrt();
        (mean, se)
    }

    /// Pathwise AD valuation with a frozen exercise policy.  The stopping
    /// decision is deliberately made on primal values (the envelope theorem),
    /// while the discounted payoff remains in the generic AD domain.
    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::needless_range_loop)]
    pub fn price_with_policy_ad<T: Real>(
        &self,
        policy: &ExerciseBoundaryPolicy,
        spot: T,
        strike: T,
        rate: T,
        div: T,
        vol: T,
        expiry: T,
    ) -> T {
        let n = self.n_steps.max(1);
        let dt = expiry / T::from_f64(n as f64);
        let rng = PhiloxRng::new(self.seed + 1_000_000);
        let mut total = T::zero();
        for p in 0..self.n_valuation_paths.max(1) {
            let mut s = spot;
            let mut payoff = T::zero();
            for j in 1..=n {
                let z = rng.normal4(p as u64, j as u32)[0];
                s *= ((rate - div - vol * vol * 0.5) * dt + vol * dt.sqrt() * z).exp();
                let exercise = (strike - s).max_val(T::zero());
                let cont = if j < n && !policy.coeffs[j].is_empty() {
                    let x = s.val() / strike.val() - 1.0;
                    policy.coeffs[j]
                        .iter()
                        .zip(basis(x, self.basis_degree))
                        .map(|(a, b)| a * b)
                        .sum::<f64>()
                } else {
                    f64::NEG_INFINITY
                };
                if j == n || (exercise.val() > cont && exercise.val() > 0.0) {
                    payoff = exercise * (-rate * dt * j as f64).exp();
                    break;
                }
            }
            total += payoff;
        }
        total / self.n_valuation_paths.max(1) as f64
    }
    pub fn price_american(
        &self,
        strike: f64,
        spot: f64,
        vol: f64,
        rate: f64,
        expiry: f64,
        is_call: bool,
    ) -> (f64, f64) {
        let p = self.calibrate_policy(strike, spot, vol, rate, expiry, is_call);
        self.price_with_policy(&p, strike, spot, vol, rate, expiry, is_call)
    }
}
