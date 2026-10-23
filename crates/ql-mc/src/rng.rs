use core::simd::Simd;
use ql_math::inv_norm_cdf;
#[derive(Clone, Copy, Debug)]
pub struct PhiloxRng {
    pub seed: u64,
}
impl PhiloxRng {
    pub fn new(seed: u64) -> Self {
        Self { seed }
    }
    fn philox(mut c: [u32; 4], mut k: [u32; 2]) -> [u32; 4] {
        for _ in 0..10 {
            let p0 = (c[0] as u64).wrapping_mul(0xD251_1F53);
            let p1 = (c[2] as u64).wrapping_mul(0xCD9E_8D57);
            c = [
                (p1 >> 32) as u32 ^ c[1] ^ k[0],
                p1 as u32,
                (p0 >> 32) as u32 ^ c[3] ^ k[1],
                p0 as u32,
            ];
            k[0] = k[0].wrapping_add(0x9E37_79B9);
            k[1] = k[1].wrapping_add(0xBB67_AE85);
        }
        c
    }
    pub fn uniform4(&self, path: u64, step: u32) -> [f64; 4] {
        let x = Self::philox(
            [path as u32, (path >> 32) as u32, step, 0],
            [self.seed as u32, (self.seed >> 32) as u32],
        );
        x.map(|v| (v as f64 + 0.5) * (1. / 4294967296.))
    }
    pub fn normal4(&self, path: u64, step: u32) -> [f64; 4] {
        let u = self.uniform4(path, step);
        // Two independent Box--Muller transforms consume all four Philox
        // uniforms.  This gives every returned component the same exact
        // Gaussian law (rather than mixing a pair with inverse-CDF draws).
        let r0 = (-2.0 * u[0].ln()).sqrt();
        let t0 = std::f64::consts::TAU * u[1];
        let r1 = (-2.0 * u[2].ln()).sqrt();
        let t1 = std::f64::consts::TAU * u[3];
        [r0 * t0.cos(), r0 * t0.sin(), r1 * t1.cos(), r1 * t1.sin()]
    }
    pub fn normal_simd8(&self, path: u64, step: u32) -> Simd<f64, 8> {
        use core::simd::Select;
        use core::simd::prelude::{SimdFloat, SimdPartialEq, SimdPartialOrd};

        let a = Self::philox(
            [path as u32, (path >> 32) as u32, step, 0],
            [self.seed as u32, (self.seed >> 32) as u32],
        );
        let b = Self::philox(
            [(path + 1) as u32, ((path + 1) >> 32) as u32, step, 0],
            [self.seed as u32, (self.seed >> 32) as u32],
        );

        const U32_TO_F64_MAGIC_BITS: u64 = 0x4330_0000_0000_0000; // 2^52
        const U32_TO_F64_MAGIC: f64 = 4_503_599_627_370_496.0; // 2^52 - 0.5
        let scale = Simd::<f64, 4>::splat(1.0 / 4294967296.0);

        let rad_u64 =
            Simd::<u64, 4>::from_array([a[0] as u64, a[2] as u64, b[0] as u64, b[2] as u64]);
        let ang_u64 =
            Simd::<u64, 4>::from_array([a[1] as u64, a[3] as u64, b[1] as u64, b[3] as u64]);
        let u_rad = (Simd::<f64, 4>::from_bits(rad_u64 | Simd::splat(U32_TO_F64_MAGIC_BITS))
            - Simd::splat(U32_TO_F64_MAGIC - 0.5))
            * scale;
        let u_ang = (Simd::<f64, 4>::from_bits(ang_u64 | Simd::splat(U32_TO_F64_MAGIC_BITS))
            - Simd::splat(U32_TO_F64_MAGIC - 0.5))
            * scale;

        // Pure-SIMD Box-Muller: r = sqrt(-2 * ln(u_rad)) without i64<->f64 casts
        let bits = u_rad.to_bits();
        let exp_u64 = (bits >> Simd::splat(52)) & Simd::splat(0x7ff);
        let mut exp_f64 = Simd::<f64, 4>::from_bits(exp_u64 | Simd::splat(U32_TO_F64_MAGIC_BITS))
            - Simd::splat(U32_TO_F64_MAGIC + 1023.0);
        let m_bits =
            (bits & Simd::splat(0x000f_ffff_ffff_ffff)) | Simd::splat(0x3ff0_0000_0000_0000);
        let mut m = Simd::<f64, 4>::from_bits(m_bits);
        let gt_sqrt2 = m.simd_gt(Simd::splat(std::f64::consts::SQRT_2));
        m = gt_sqrt2.select(m * Simd::splat(0.5), m);
        exp_f64 = gt_sqrt2.select(exp_f64 + Simd::splat(1.0), exp_f64);

        let z = (m - Simd::splat(1.0)) / (m + Simd::splat(1.0));
        let z2 = z * z;
        let z4 = z2 * z2;
        let p01 = Simd::splat(1.0) + z2 * Simd::splat(1.0 / 3.0);
        let p23 = Simd::splat(1.0 / 5.0) + z2 * Simd::splat(1.0 / 7.0);
        let p45 = Simd::splat(1.0 / 9.0) + z2 * Simd::splat(1.0 / 11.0);
        let p67 = Simd::splat(1.0 / 13.0) + z2 * Simd::splat(1.0 / 15.0);
        let ln_p = (p01 + z4 * p23) + (z4 * z4) * (p45 + z4 * p67);
        let ln_u = exp_f64 * Simd::splat(std::f64::consts::LN_2) + Simd::splat(2.0) * z * ln_p;
        let rad = std::simd::StdFloat::sqrt(Simd::splat(-2.0) * ln_u);

        // Pure-SIMD sin(2*pi*u_ang) and cos(2*pi*u_ang) with Estrin evaluation
        const MAGIC: f64 = 6_755_399_441_055_744.0;
        let u_q = u_ang * Simd::splat(4.0) + Simd::splat(MAGIC);
        let q_f64 = u_q - Simd::splat(MAGIC);
        let r = (u_ang - q_f64 * Simd::splat(0.25)) * Simd::splat(std::f64::consts::TAU);
        let r2 = r * r;
        let r4 = r2 * r2;

        let s01 = Simd::splat(1.0) + r2 * Simd::splat(-1.0 / 6.0);
        let s23 = Simd::splat(1.0 / 120.0) + r2 * Simd::splat(-1.0 / 5_040.0);
        let s45 = Simd::splat(1.0 / 362_880.0) + r2 * Simd::splat(-1.0 / 39_916_800.0);
        let sin_r = r * ((s01 + r4 * s23) + (r4 * r4) * s45);

        let c01 = Simd::splat(1.0) + r2 * Simd::splat(-0.5);
        let c23 = Simd::splat(1.0 / 24.0) + r2 * Simd::splat(-1.0 / 720.0);
        let c45 = Simd::splat(1.0 / 40_320.0) + r2 * Simd::splat(-1.0 / 3_628_800.0);
        let cos_r = (c01 + r4 * c23) + (r4 * r4) * c45;

        let qi = u_q.to_bits() & Simd::splat(3);
        let swap_mask = (qi & Simd::splat(1)).simd_ne(Simd::splat(0));
        let flip_mask = (qi & Simd::splat(2)).simd_ne(Simd::splat(0));
        let mut cos_t = swap_mask.select(-sin_r, cos_r);
        let mut sin_t = swap_mask.select(cos_r, sin_r);
        cos_t = flip_mask.select(-cos_t, cos_t);
        sin_t = flip_mask.select(-sin_t, sin_t);

        let n_cos = (rad * cos_t).to_array();
        let n_sin = (rad * sin_t).to_array();
        Simd::from_array([
            n_cos[0], n_sin[0], n_cos[1], n_sin[1], n_cos[2], n_sin[2], n_cos[3], n_sin[3],
        ])
    }
}
#[derive(Clone, Debug)]
pub struct SobolGenerator {
    pub dimension: usize,
    pub index: u64,
    pub state: Vec<u32>,
    pub dir_nums: Vec<[u32; 32]>,
}
impl SobolGenerator {
    #[allow(clippy::needless_range_loop)]
    pub fn new(dimension: usize) -> Self {
        assert!(
            dimension <= 64,
            "SobolGenerator supports at most 64 dimensions"
        );
        let d = dimension;
        // Joe--Kuo new-joe-kuo-6.21201 primitive polynomials and initial
        // direction numerators for dimensions 2..64.  Dimension one is the
        // Van der Corput sequence below.  Keeping the actual published seeds
        // (rather than synthesizing arbitrary m_k values) makes the generator
        // reproducible with standard Sobol implementations.
        const PARAMS: &[(usize, u32, &[u32])] = &[
            (1, 0, &[1]),
            (2, 1, &[1, 3]),
            (3, 1, &[1, 3, 1]),
            (3, 2, &[1, 1, 1]),
            (4, 1, &[1, 1, 3, 3]),
            (4, 4, &[1, 3, 5, 13]),
            (5, 2, &[1, 1, 5, 5, 17]),
            (5, 4, &[1, 1, 5, 5, 5]),
            (5, 7, &[1, 1, 7, 11, 19]),
            (5, 11, &[1, 1, 5, 1, 1]),
            (5, 13, &[1, 1, 1, 3, 11]),
            (5, 14, &[1, 3, 5, 5, 31]),
            (6, 1, &[1, 3, 3, 9, 7, 49]),
            (6, 13, &[1, 1, 1, 15, 21, 21]),
            (6, 16, &[1, 3, 1, 13, 27, 49]),
            (6, 19, &[1, 1, 1, 15, 7, 5]),
            (6, 22, &[1, 3, 1, 15, 13, 25]),
            (6, 25, &[1, 1, 5, 5, 19, 61]),
            (7, 1, &[1, 3, 7, 11, 23, 15, 103]),
            (7, 4, &[1, 3, 7, 13, 13, 15, 69]),
            (7, 7, &[1, 1, 3, 13, 7, 35, 63]),
            (7, 8, &[1, 3, 5, 9, 1, 25, 53]),
            (7, 14, &[1, 3, 1, 13, 9, 35, 107]),
            (7, 19, &[1, 3, 1, 5, 27, 61, 31]),
            (7, 21, &[1, 1, 5, 11, 19, 41, 61]),
            (7, 28, &[1, 3, 5, 3, 3, 13, 69]),
            (7, 31, &[1, 1, 7, 13, 1, 19, 1]),
            (7, 32, &[1, 3, 7, 5, 13, 19, 59]),
            (7, 37, &[1, 1, 3, 9, 25, 29, 41]),
            (7, 41, &[1, 3, 5, 13, 23, 1, 55]),
            (7, 42, &[1, 3, 7, 3, 13, 59, 17]),
            (7, 50, &[1, 3, 1, 3, 5, 53, 69]),
            (7, 55, &[1, 1, 5, 5, 23, 33, 13]),
            (7, 56, &[1, 1, 7, 7, 1, 61, 123]),
            (7, 59, &[1, 1, 7, 9, 13, 61, 49]),
            (7, 62, &[1, 3, 3, 5, 3, 55, 33]),
            (8, 14, &[1, 3, 1, 15, 31, 13, 49, 245]),
            (8, 21, &[1, 3, 5, 15, 31, 59, 63, 97]),
            (8, 22, &[1, 3, 1, 11, 11, 11, 77, 249]),
            (8, 38, &[1, 3, 1, 11, 27, 43, 71, 9]),
            (8, 47, &[1, 1, 7, 15, 21, 11, 81, 45]),
            (8, 49, &[1, 3, 7, 3, 25, 31, 65, 79]),
            (8, 50, &[1, 3, 1, 1, 19, 11, 3, 205]),
            (8, 52, &[1, 1, 5, 9, 19, 21, 29, 157]),
            (8, 56, &[1, 3, 7, 11, 1, 33, 89, 185]),
            (8, 67, &[1, 3, 3, 3, 15, 9, 79, 71]),
            (8, 70, &[1, 3, 7, 11, 15, 39, 119, 27]),
            (8, 84, &[1, 1, 3, 1, 11, 31, 97, 225]),
            (8, 97, &[1, 1, 1, 3, 23, 43, 57, 177]),
            (8, 103, &[1, 3, 7, 7, 17, 17, 37, 71]),
            (8, 115, &[1, 3, 1, 5, 27, 63, 123, 213]),
            (8, 122, &[1, 1, 3, 5, 11, 43, 53, 133]),
            (9, 8, &[1, 3, 5, 5, 29, 17, 47, 173, 479]),
            (9, 13, &[1, 3, 3, 11, 3, 1, 109, 9, 69]),
            (9, 16, &[1, 1, 1, 5, 17, 39, 23, 5, 343]),
            (9, 22, &[1, 3, 1, 5, 25, 15, 31, 103, 499]),
            (9, 25, &[1, 1, 1, 11, 11, 17, 63, 105, 183]),
            (9, 44, &[1, 1, 5, 11, 9, 29, 97, 231, 363]),
            (9, 47, &[1, 1, 5, 15, 19, 45, 41, 7, 383]),
            (9, 52, &[1, 3, 7, 7, 31, 19, 83, 137, 221]),
            (9, 55, &[1, 1, 1, 3, 23, 15, 111, 223, 83]),
            (9, 59, &[1, 1, 5, 13, 31, 15, 55, 25, 161]),
            (9, 62, &[1, 1, 3, 13, 25, 47, 39, 87, 257]),
        ];
        let mut dirs = vec![[0u32; 32]; d];
        if d > 0 {
            for b in 0..32 {
                dirs[0][b] = 1u32 << (31 - b);
            }
        }
        for dim in 1..d {
            let (deg, poly, initial) = PARAMS[dim - 1];
            let mut v = [0u32; 32];
            for k in 1..=32 {
                if k <= deg {
                    v[k - 1] = initial[k - 1] << (32 - k);
                } else {
                    let mut x = v[k - 1 - deg] ^ (v[k - 1 - deg] >> deg);
                    for j in 1..deg {
                        if ((poly >> (deg - 1 - j)) & 1) != 0 {
                            x ^= v[k - 1 - j];
                        }
                    }
                    v[k - 1] = x;
                }
            }
            dirs[dim] = v;
        }
        Self {
            dimension: d,
            index: 0,
            state: vec![0; d],
            dir_nums: dirs,
        }
    }
    pub fn next_point(&mut self, out: &mut [f64]) {
        self.index += 1;
        let c = (self.index.trailing_zeros().min(31)) as usize;
        for d in 0..self.dimension.min(out.len()) {
            self.state[d] ^= self.dir_nums[d][c];
            out[d] = (self.state[d] as f64 + 0.5) * (1. / 4294967296.);
        }
    }
    pub fn next_normal_point(&mut self, out: &mut [f64]) {
        self.next_point(out);
        let n = out.len().min(self.dimension);
        for x in &mut out[..n] {
            *x = inv_norm_cdf((*x).clamp(1e-15, 1. - 1e-15));
        }
    }

    /// Stateless random-access Sobol integer point.  The Gray-code XOR makes
    /// this exactly agree with the sequential generator at the same index.
    pub fn point_u32_at(&self, index: u32, out: &mut [u32]) {
        let dim = self.dimension.min(out.len());
        out[..dim].fill(0);
        let mut bits = index ^ (index >> 1);
        while bits != 0 {
            let bit = bits.trailing_zeros() as usize;
            for d in 0..dim {
                out[d] ^= self.dir_nums[d][bit];
            }
            bits &= bits - 1;
        }
    }
}

#[inline]
pub fn laine_karras_permutation(mut x: u32, seed: u32) -> u32 {
    x = x.wrapping_add(seed);
    x ^= x.wrapping_mul(0x6c50_b47c);
    x ^= x.wrapping_mul(0xb82f_1e52);
    x ^= x.wrapping_mul(0xc7af_e638);
    x ^= x.wrapping_mul(0x8d22_f6e6);
    x
}

#[inline]
pub fn nested_uniform_scramble(mut x: u32, seed: u32) -> u32 {
    x = x.reverse_bits();
    x = laine_karras_permutation(x, seed);
    x.reverse_bits()
}

#[inline]
pub fn hash_combine_u32(a: u32, b: u32) -> u32 {
    let mut h = a ^ b
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(a << 6)
        .wrapping_add(a >> 2);
    h ^= h >> 16;
    h = h.wrapping_mul(0x85eb_ca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2_ae35);
    h ^ (h >> 16)
}

#[derive(Clone, Debug)]
pub struct BurleyOwenSobol {
    sobol: SobolGenerator,
    seed: u32,
    index: u32,
    dim_seeds: [u32; 64],
}
impl BurleyOwenSobol {
    pub fn new(dimension: usize, seed: u32) -> Self {
        let mut dim_seeds = [0u32; 64];
        for (d, slot) in dim_seeds.iter_mut().enumerate().take(dimension.min(64)) {
            *slot = hash_combine_u32(seed, d as u32 + 1);
        }
        Self {
            sobol: SobolGenerator::new(dimension),
            seed,
            index: 0,
            dim_seeds,
        }
    }
    pub fn dimension(&self) -> usize {
        self.sobol.dimension
    }
    pub fn point_at(&self, index: u32, out: &mut [f64]) {
        let mut raw = [0u32; 64];
        let shuffled = nested_uniform_scramble(index, self.seed);
        self.sobol
            .point_u32_at(shuffled, &mut raw[..self.sobol.dimension]);
        for d in 0..self.sobol.dimension.min(out.len()) {
            let x = nested_uniform_scramble(raw[d], self.dim_seeds[d]);
            out[d] = (x as f64 + 0.5) * (1.0 / 4_294_967_296.0);
        }
    }
    pub fn normal_point_at(&self, index: u32, out: &mut [f64]) {
        self.point_at(index, out);
        for x in out.iter_mut() {
            *x = ql_math::inv_norm_cdf((*x).clamp(1e-15, 1.0 - 1e-15));
        }
    }
    pub fn next_point(&mut self, out: &mut [f64]) {
        self.point_at(self.index, out);
        self.index = self.index.wrapping_add(1);
    }
    pub fn next_normal_point(&mut self, out: &mut [f64]) {
        self.normal_point_at(self.index, out);
        self.index = self.index.wrapping_add(1);
    }
}

pub fn rqmc_estimate_with_stderr<F>(
    dimension: usize,
    points_per_scramble: u32,
    scrambles: u32,
    base_seed: u32,
    integrand: F,
) -> (f64, f64)
where
    F: Fn(&[f64]) -> f64,
{
    let m = scrambles.max(2);
    let n = points_per_scramble.max(1);
    let mut means = Vec::with_capacity(m as usize);
    let mut point = vec![0.0; dimension];
    for s in 0..m {
        let sobol = BurleyOwenSobol::new(dimension, hash_combine_u32(base_seed, s));
        let mut sum = 0.0;
        for i in 0..n {
            sobol.normal_point_at(i, &mut point);
            sum += integrand(&point);
        }
        means.push(sum / n as f64);
    }
    let mean = means.iter().sum::<f64>() / m as f64;
    let var = means.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (m.saturating_sub(1) as f64);
    (mean, (var / m as f64).sqrt())
}
#[derive(Clone, Debug)]
pub struct BrownianBridge {
    pub times: Vec<f64>,
    steps: Vec<(Option<usize>, usize, usize, f64, f64, f64)>,
}
impl BrownianBridge {
    pub fn new(times: Vec<f64>) -> Self {
        assert!(
            times.iter().all(|t| t.is_finite() && *t > 0.0)
                && times.windows(2).all(|w| w[1] > w[0]),
            "Brownian-bridge times must be strictly increasing and positive"
        );
        let mut steps = Vec::new();
        let n = times.len();
        let mut q = std::collections::VecDeque::new();
        if n > 1 {
            q.push_back((None, n - 1));
            while let Some((l, r)) = q.pop_front() {
                let start = l.map_or(0, |x| x + 1);
                if start >= r {
                    continue;
                }
                let m = (start + r - 1) / 2;
                let tl = l.map_or(0., |x| times[x]);
                let tm = times[m];
                let tr = times[r];
                let den = tr - tl;
                steps.push((
                    l,
                    r,
                    m,
                    (tr - tm) / den,
                    (tm - tl) / den,
                    ((tm - tl) * (tr - tm) / den).sqrt(),
                ));
                q.push_back((l, m));
                q.push_back((Some(m), r));
            }
        }
        Self { times, steps }
    }
    pub fn transform(&self, normals: &[f64], out: &mut [f64]) {
        let n = self.times.len();
        if n == 0 {
            return;
        }
        assert!(
            normals.len() >= n && out.len() >= n,
            "BrownianBridge::transform requires one normal and output slot per time"
        );
        let mut w_stack = [0.0f64; 64];
        let mut w_heap = Vec::new();
        let w: &mut [f64] = if n <= 64 {
            &mut w_stack[..n]
        } else {
            w_heap.resize(n, 0.0);
            &mut w_heap
        };
        w[n - 1] = self.times[n - 1].sqrt() * normals[0];
        for (k, &(l, r, m, wl, wr, sig)) in self.steps.iter().enumerate() {
            if k + 1 >= normals.len() {
                break;
            }
            let left = l.map_or(0., |i| w[i]);
            w[m] = wl * left + wr * w[r] + sig * normals[k + 1];
        }
        let mut prev = 0.;
        for i in 0..n {
            out[i] = w[i] - prev;
            prev = w[i];
        }
    }
}
