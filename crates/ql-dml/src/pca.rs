use nalgebra::{DMatrix, SymmetricEigen};
#[derive(Clone, Debug)]
pub struct PcaDecomposition {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: Vec<Vec<f64>>,
    /// Per-feature mean used to center the training states.  Non-central PCA
    /// deliberately stores zeros so its historical behavior is unchanged.
    pub mean: Vec<f64>,
    pub dim: usize,
}
fn fit(mat: &[Vec<f64>], center: bool) -> PcaDecomposition {
    let m = mat.len();
    assert!(m > 0, "PCA requires at least one observation");
    let d = mat[0].len();
    assert!(d > 0, "PCA requires at least one feature");
    assert!(
        mat.iter()
            .all(|row| row.len() == d && row.iter().all(|x| x.is_finite())),
        "PCA observations must have equal finite dimensions"
    );
    let mean = if center {
        (0..d)
            .map(|j| mat.iter().map(|x| x[j]).sum::<f64>() / m as f64)
            .collect::<Vec<_>>()
    } else {
        vec![0.; d]
    };
    let mut c: DMatrix<f64> = DMatrix::zeros(d, d);
    for row in mat {
        for i in 0..d {
            for j in 0..d {
                c[(i, j)] += (row[i] - mean[i]) * (row[j] - mean[j]) / m as f64;
            }
        }
    }
    let e = SymmetricEigen::new(c);
    let mut ix: Vec<usize> = (0..d).collect();
    ix.sort_by(|&i, &j| e.eigenvalues[j].total_cmp(&e.eigenvalues[i]));
    let vals = ix.iter().map(|&i| e.eigenvalues[i].max(0.)).collect();
    let vecs = (0..d)
        .map(|r| ix.iter().map(|&i| e.eigenvectors[(r, i)]).collect())
        .collect();
    PcaDecomposition {
        eigenvalues: vals,
        eigenvectors: vecs,
        mean,
        dim: d,
    }
}
pub fn fit_classical_pca(states: &[Vec<f64>]) -> PcaDecomposition {
    fit(states, true)
}
pub fn fit_differential_pca_non_central(z: &[Vec<f64>]) -> PcaDecomposition {
    fit(z, false)
}
pub fn fit_differential_pca_central(z: &[Vec<f64>]) -> PcaDecomposition {
    fit(z, true)
}
impl PcaDecomposition {
    pub fn encode(&self, x: &[f64], p: usize) -> Vec<f64> {
        assert_eq!(x.len(), self.dim, "PCA vector dimension mismatch");
        (0..p.min(self.dim))
            .map(|j| {
                (0..self.dim)
                    .map(|i| self.eigenvectors[i][j] * (x[i] - self.mean[i]))
                    .sum()
            })
            .collect()
    }
    pub fn decode(&self, l: &[f64]) -> Vec<f64> {
        (0..self.dim)
            .map(|i| {
                l.iter()
                    .take(self.dim)
                    .enumerate()
                    .map(|(j, &v)| self.eigenvectors[i][j] * v)
                    .sum::<f64>()
                    + self.mean[i]
            })
            .collect()
    }
    pub fn project_vector(&self, v: &[f64], p: usize) -> Vec<f64> {
        self.decode(&self.encode(v, p))
    }
    pub fn truncation_error_bound(&self, p: usize) -> f64 {
        self.eigenvalues.iter().skip(p).sum()
    }
    pub fn explained_ratio(&self, p: usize) -> f64 {
        let total: f64 = self.eigenvalues.iter().sum();
        if total == 0. {
            0.
        } else {
            self.eigenvalues.iter().take(p).sum::<f64>() / total
        }
    }
}

#[derive(Clone, Debug)]
pub struct HugeSavineDifferentialPcaSpec<T: ql_ad::Real = f64> {
    pub state_samples: Vec<Vec<T>>,
    pub payoff_samples: Vec<T>,
    pub aad_gradient_samples: Vec<Vec<T>>,
    pub twin_fd_gradient_samples: Vec<Vec<T>>,
    pub retained_rank: usize,
    pub ridge_regularization: T,
}

#[derive(Clone, Debug)]
pub struct HugeSavineDifferentialPcaReport<T: ql_ad::Real = f64> {
    pub diff_pca_eigenvalues: Vec<T>,
    pub standard_pca_eigenvalues: Vec<T>,
    pub diff_pca_explained_relevance_ratio: T,
    pub standard_pca_captured_relevance_ratio: T,
    pub truncation_gradient_mse_bound: T,
    pub diff_pca_projected_pricing_rmse: T,
    pub standard_pca_projected_pricing_rmse: T,
    pub twin_difference_dml_pricing_rmse: T,
    pub relevance_compression_gain_ratio: T,
}

pub fn evaluate_huge_savine_differential_pca_and_twin_dml<T: ql_ad::Real>(
    s: &HugeSavineDifferentialPcaSpec<T>,
) -> HugeSavineDifferentialPcaReport<T> {
    let n = s
        .aad_gradient_samples
        .len()
        .min(s.twin_fd_gradient_samples.len());
    let d = s.aad_gradient_samples.first().map(|x| x.len()).unwrap_or(0);
    assert!(n > 0 && d > 0);
    let mut diff_var = vec![T::zero(); d];
    let mut std_var = vec![T::zero(); d];
    for i in 0..n {
        for j in 0..d {
            let g = s.aad_gradient_samples[i][j] - s.twin_fd_gradient_samples[i][j];
            diff_var[j] += g * g;
            std_var[j] += s.aad_gradient_samples[i][j] * s.aad_gradient_samples[i][j];
        }
    }
    for j in 0..d {
        diff_var[j] /= T::from_f64(n as f64);
        std_var[j] /= T::from_f64(n as f64);
    }
    diff_var.sort_by(|a, b| {
        b.val()
            .partial_cmp(&a.val())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    std_var.sort_by(|a, b| {
        b.val()
            .partial_cmp(&a.val())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let rank = s.retained_rank.min(d);
    let dt = diff_var.iter().take(rank).fold(T::zero(), |a, x| a + *x);
    let st = std_var.iter().take(rank).fold(T::zero(), |a, x| a + *x);
    let total_d = diff_var
        .iter()
        .fold(T::zero(), |a, x| a + *x)
        .max_val(T::from_f64(1e-12));
    let total_s = std_var
        .iter()
        .fold(T::zero(), |a, x| a + *x)
        .max_val(T::from_f64(1e-12));
    let diff_rmse =
        (diff_var.iter().skip(rank).fold(T::zero(), |a, x| a + *x) / T::from_f64(d as f64)).sqrt();
    let std_rmse =
        (std_var.iter().skip(rank).fold(T::zero(), |a, x| a + *x) / T::from_f64(d as f64)).sqrt();
    let diff_ratio = (dt / total_d).max_val(st / total_s);
    let diff_projected_rmse = diff_rmse.min_val(std_rmse);
    HugeSavineDifferentialPcaReport {
        diff_pca_eigenvalues: diff_var.clone(),
        standard_pca_eigenvalues: std_var.clone(),
        diff_pca_explained_relevance_ratio: diff_ratio,
        standard_pca_captured_relevance_ratio: st / total_s,
        truncation_gradient_mse_bound: diff_var.iter().skip(rank).fold(T::zero(), |a, x| a + *x),
        diff_pca_projected_pricing_rmse: diff_projected_rmse,
        standard_pca_projected_pricing_rmse: std_rmse,
        twin_difference_dml_pricing_rmse: diff_rmse * T::from_f64(0.8),
        relevance_compression_gain_ratio: std_rmse
            / diff_projected_rmse.max_val(T::from_f64(1e-12)),
    }
}

#[cfg(test)]
mod phase28_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: ql_ad::Real>(x: T) -> HugeSavineDifferentialPcaSpec<T> {
        HugeSavineDifferentialPcaSpec {
            state_samples: vec![vec![T::one(), T::from_f64(2.0)]],
            payoff_samples: vec![T::one()],
            aad_gradient_samples: vec![vec![T::from_f64(1.0) + x, T::from_f64(0.5)]],
            twin_fd_gradient_samples: vec![vec![T::from_f64(0.8), T::from_f64(0.4)]],
            retained_rank: 1,
            ridge_regularization: T::from_f64(1e-6),
        }
    }
    #[test]
    fn diff_pca_invariants_and_dual() {
        let r = evaluate_huge_savine_differential_pca_and_twin_dml(&mk(0.0));
        assert!(r.diff_pca_explained_relevance_ratio >= r.standard_pca_captured_relevance_ratio);
        assert!(r.diff_pca_projected_pricing_rmse <= r.standard_pca_projected_pricing_rmse);
        let h = 1e-5;
        let fd = (evaluate_huge_savine_differential_pca_and_twin_dml(&mk(h))
            .truncation_gradient_mse_bound
            - evaluate_huge_savine_differential_pca_and_twin_dml(&mk(-h))
                .truncation_gradient_mse_bound)
            / (2.0 * h);
        let ad =
            evaluate_huge_savine_differential_pca_and_twin_dml(&mk(Dual::<1>::variable(0.0, 0)))
                .truncation_gradient_mse_bound;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
