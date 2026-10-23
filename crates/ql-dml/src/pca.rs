use nalgebra::{DMatrix, SymmetricEigen};
#[derive(Clone, Debug)]
pub struct PcaDecomposition {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: Vec<Vec<f64>>,
    pub dim: usize,
}
fn fit(mat: &[Vec<f64>], center: bool) -> PcaDecomposition {
    let m = mat.len();
    assert!(m > 0);
    let d = mat[0].len();
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
        (0..p.min(self.dim))
            .map(|j| (0..self.dim).map(|i| self.eigenvectors[i][j] * x[i]).sum())
            .collect()
    }
    pub fn decode(&self, l: &[f64]) -> Vec<f64> {
        (0..self.dim)
            .map(|i| {
                l.iter()
                    .enumerate()
                    .map(|(j, &v)| self.eigenvectors[i][j] * v)
                    .sum()
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
