//! Counter-based Philox4x32-10, shared by model and path engines.
use core::simd::Simd;

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
