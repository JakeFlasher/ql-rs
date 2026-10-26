#![allow(clippy::needless_range_loop)]
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

fn sym_eigen_jacobi<T: ql_ad::Real>(mat: &[Vec<T>], sweeps: usize) -> (Vec<T>, Vec<Vec<T>>) {
    let n = mat.len();
    let mut a = mat.to_vec();
    let mut v = vec![vec![T::zero(); n]; n];
    for i in 0..n {
        v[i][i] = T::one();
    }
    for _ in 0..sweeps.max(2) {
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p][q];
                let diff = a[q][q] - a[p][p];
                let disc = (diff * diff + T::from_f64(4.0) * apq * apq + T::from_f64(1e-20)).sqrt();
                let c =
                    (T::one() + diff / disc).max_val(T::zero()).sqrt() / T::from_f64(2.0).sqrt();
                let sign = if apq.val() >= 0.0 {
                    T::one()
                } else {
                    -T::one()
                };
                let sn = sign * (T::one() - c * c).max_val(T::zero()).sqrt();
                for k in 0..n {
                    let akp = a[k][p];
                    let akq = a[k][q];
                    a[k][p] = c * akp - sn * akq;
                    a[k][q] = sn * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p][k];
                    let aqk = a[q][k];
                    a[p][k] = c * apk - sn * aqk;
                    a[q][k] = sn * apk + c * aqk;
                    let vkp = v[k][p];
                    let vkq = v[k][q];
                    v[k][p] = c * vkp - sn * vkq;
                    v[k][q] = sn * vkp + c * vkq;
                }
            }
        }
    }
    let mut ix: Vec<usize> = (0..n).collect();
    ix.sort_by(|&i, &j| {
        a[j][j]
            .val()
            .partial_cmp(&a[i][i].val())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let vals = ix
        .iter()
        .map(|&i| a[i][i].max_val(T::zero()))
        .collect::<Vec<_>>();
    let vecs = (0..n)
        .map(|r| ix.iter().map(|&i| v[r][i]).collect())
        .collect();
    (vals, vecs)
}

pub fn evaluate_huge_savine_differential_pca_and_twin_dml<T: ql_ad::Real>(
    s: &HugeSavineDifferentialPcaSpec<T>,
) -> HugeSavineDifferentialPcaReport<T> {
    let n = s
        .state_samples
        .len()
        .min(s.payoff_samples.len())
        .min(s.aad_gradient_samples.len())
        .min(s.twin_fd_gradient_samples.len());
    assert!(n > 0);
    let d = s.state_samples[0].len();
    assert!(d > 0);
    let mut cg = vec![vec![T::zero(); d]; d];
    let mut sx = vec![vec![T::zero(); d]; d];
    for i in 0..n {
        for p in 0..d {
            for q in 0..d {
                cg[p][q] += s.aad_gradient_samples[i][p] * s.aad_gradient_samples[i][q]
                    / T::from_f64(n as f64);
                sx[p][q] += s.state_samples[i][p] * s.state_samples[i][q] / T::from_f64(n as f64);
            }
        }
    }
    let ridge = s.ridge_regularization.max_val(T::zero()) * T::from_f64(1e-6);
    for i in 0..d {
        cg[i][i] += ridge;
        sx[i][i] += ridge;
    }
    let (eval_d, vec_d) = sym_eigen_jacobi(&cg, 8 * d * d);
    let (eval_s, vec_s) = sym_eigen_jacobi(&sx, 8 * d * d);
    let rank = s.retained_rank.min(d);
    let total_d = eval_d
        .iter()
        .copied()
        .fold(T::zero(), |a, x| a + x)
        .max_val(T::from_f64(1e-12));
    let mut diff_sum = T::zero();
    for x in eval_d.iter().take(rank) {
        diff_sum += *x;
    }
    let mut std_cap = T::zero();
    for j in 0..rank {
        let mut q = T::zero();
        for p in 0..d {
            for r in 0..d {
                q += vec_s[p][j] * cg[p][r] * vec_s[r][j];
            }
        }
        std_cap += q;
    }
    let mut tr = T::zero();
    for x in eval_d.iter() {
        tr += *x;
    }
    let mut diff_mse = T::zero();
    for x in eval_d.iter().skip(rank) {
        diff_mse += *x;
    }
    let trunc_bound = diff_mse;
    let mut standard_mse = T::zero();
    let mut twin_mse = T::zero();
    let mut proj_diff_mse = T::zero();
    for i in 0..n {
        let mut zd = vec![T::zero(); rank];
        let mut zs = vec![T::zero(); rank];
        for j in 0..rank {
            for p in 0..d {
                zd[j] += vec_d[p][j] * s.state_samples[i][p];
                zs[j] += vec_s[p][j] * s.state_samples[i][p];
            }
        }
        let mut pd = T::zero();
        let mut ps = T::zero();
        for j in 0..rank {
            let mut numd = T::zero();
            let mut nums = T::zero();
            let mut den_d = T::zero();
            let mut den_s = T::zero();
            for k in 0..n {
                let mut zdk = T::zero();
                let mut zsk = T::zero();
                for p in 0..d {
                    zdk += vec_d[p][j] * s.state_samples[k][p];
                    zsk += vec_s[p][j] * s.state_samples[k][p];
                }
                numd += zdk * s.payoff_samples[k];
                nums += zsk * s.payoff_samples[k];
                den_d += zdk * zdk + s.ridge_regularization;
                den_s += zsk * zsk + s.ridge_regularization;
            }
            pd += zd[j] * numd / den_d.max_val(T::from_f64(1e-12));
            ps += zs[j] * nums / den_s.max_val(T::from_f64(1e-12));
        }
        let mut twin = T::zero();
        for p in 0..d {
            twin += s.state_samples[i][p]
                * (s.aad_gradient_samples[i][p] + s.twin_fd_gradient_samples[i][p])
                * T::from_f64(0.5);
        }
        let epd = pd - s.payoff_samples[i];
        let eps = ps - s.payoff_samples[i];
        let et = twin - s.payoff_samples[i];
        proj_diff_mse += epd * epd;
        standard_mse += eps * eps;
        twin_mse += et * et;
    }
    let diff_rmse = (proj_diff_mse / T::from_f64(n.max(1) as f64)).sqrt();
    let std_rmse = (standard_mse / T::from_f64(n.max(1) as f64)).sqrt();
    let twin_rmse = (twin_mse / T::from_f64(n.max(1) as f64)).sqrt();
    HugeSavineDifferentialPcaReport {
        diff_pca_eigenvalues: eval_d,
        standard_pca_eigenvalues: eval_s,
        diff_pca_explained_relevance_ratio: diff_sum / total_d,
        standard_pca_captured_relevance_ratio: std_cap / tr.max_val(T::from_f64(1e-12)),
        truncation_gradient_mse_bound: trunc_bound,
        diff_pca_projected_pricing_rmse: diff_rmse,
        standard_pca_projected_pricing_rmse: std_rmse,
        twin_difference_dml_pricing_rmse: twin_rmse,
        relevance_compression_gain_ratio: std_rmse / diff_rmse.max_val(T::from_f64(1e-12)),
    }
}

#[cfg(test)]
mod domain_suite_40_tests {
    use super::*;
    use ql_ad::Dual;
    fn mk<T: ql_ad::Real>(x: T) -> HugeSavineDifferentialPcaSpec<T> {
        HugeSavineDifferentialPcaSpec {
            state_samples: vec![
                vec![T::from_f64(0.5) + x, T::from_f64(2.0)],
                vec![T::from_f64(-0.4), T::from_f64(-1.8)],
                vec![T::from_f64(0.3), T::from_f64(1.6)],
            ],
            payoff_samples: vec![T::from_f64(1.1) + x, T::from_f64(-0.88), T::from_f64(0.66)],
            aad_gradient_samples: vec![
                vec![T::from_f64(2.2) + x, T::from_f64(0.1)],
                vec![T::from_f64(1.9), T::from_f64(-0.1)],
                vec![T::from_f64(2.0), T::from_f64(0.05)],
            ],
            twin_fd_gradient_samples: vec![
                vec![T::from_f64(2.0), T::from_f64(0.1)],
                vec![T::from_f64(1.8), T::from_f64(-0.1)],
                vec![T::from_f64(1.9), T::from_f64(0.05)],
            ],
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
