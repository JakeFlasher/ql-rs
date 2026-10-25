use crate::LevyModel;
use ql_ad::Real;
use ql_instruments::OptionType;
#[derive(Clone, Copy, Debug)]
pub struct CgmyProcess<T: Real> {
    pub sigma: T,
    pub c: T,
    pub g: T,
    pub m: T,
    pub y: T,
}
impl<T: Real> CgmyProcess<T> {
    /// Lanczos gamma approximation evaluated in the AD scalar type.  CGMY
    /// uses Gamma(2-Y)/(Y(Y-1)); keeping the approximation in `T` preserves
    /// curvature Greeks with respect to `Y` instead of lifting a finite
    /// difference constant from f64.
    fn gamma_lanczos(z: T) -> T {
        const C: [f64; 9] = [
            0.9999999999998099,
            676.5203681218851,
            -1259.1392167224028,
            771.3234287776531,
            -176.6150291621406,
            12.507343278686905,
            -0.13857109526572012,
            9.984369578019572e-6,
            1.5056327351493116e-7,
        ];
        // All production CGMY calibrations have 0 < Y < 2, hence 2-Y > 0;
        // the reflection branch is retained for a well-defined generic API.
        if z.val() < 0.5 {
            let pi = T::from_f64(std::f64::consts::PI);
            pi / ((pi * z).sin() * Self::gamma_lanczos(T::one() - z))
        } else {
            let zz = z - 1.0;
            let mut x = T::from_f64(C[0]);
            for (i, &ci) in C.iter().enumerate().skip(1) {
                x += T::from_f64(ci) / (zz + i as f64);
            }
            let t = zz + 7.5;
            T::from_f64((2.0 * std::f64::consts::PI).sqrt()) * t.pow(zz + 0.5) * (-t).exp() * x
        }
    }
    fn gamma_neg_y(&self) -> T {
        Self::gamma_lanczos(T::from_f64(2.0) - self.y) / (self.y * (self.y - 1.0))
    }
}
impl<T: Real> LevyModel<T> for CgmyProcess<T> {
    fn martingale_compensator(&self) -> T {
        -self.sigma * self.sigma * 0.5
            - self.c
                * self.gamma_neg_y()
                * ((self.m - 1.0).pow(self.y) - self.m.pow(self.y) + (self.g + 1.0).pow(self.y)
                    - self.g.pow(self.y))
    }
    fn char_exponent(&self, u: T) -> (T, T) {
        let rm = (self.m * self.m + u * u).sqrt();
        let rg = (self.g * self.g + u * u).sqrt();
        let am = (self.y * (u / self.m).atan()).cos();
        let sm = (self.y * (u / self.m).atan()).sin();
        let ag = (self.y * (u / self.g).atan()).cos();
        let sg = (self.y * (u / self.g).atan()).sin();
        let mm = rm.pow(self.y);
        let gg = rg.pow(self.y);
        let re = -self.sigma * self.sigma * u * u * 0.5
            + self.c
                * self.gamma_neg_y()
                * (mm * am - self.m.pow(self.y) + gg * ag - self.g.pow(self.y));
        let im = self.c * self.gamma_neg_y() * (-mm * sm + gg * sg);
        (re, im)
    }
    fn char_exponent_shifted(&self, u: T, beta: T) -> (T, T) {
        let mb = self.m - beta;
        let gb = self.g + beta;
        let rm = (mb * mb + u * u).sqrt();
        let rg = (gb * gb + u * u).sqrt();
        let am = (self.y * (u / mb).atan()).cos();
        let sm = (self.y * (u / mb).atan()).sin();
        let ag = (self.y * (u / gb).atan()).cos();
        let sg = (self.y * (u / gb).atan()).sin();
        let mm = rm.pow(self.y);
        let gg = rg.pow(self.y);
        let b2_u2 = beta * beta - u * u;
        let re = self.sigma * self.sigma * b2_u2 * 0.5
            + self.c
                * self.gamma_neg_y()
                * (mm * am - self.m.pow(self.y) + gg * ag - self.g.pow(self.y));
        let im =
            self.sigma * self.sigma * u * beta + self.c * self.gamma_neg_y() * (-mm * sm + gg * sg);
        (re, im)
    }
    fn cumulants(&self, r: T, q: T, t: T) -> (T, T, T) {
        let om = self.martingale_compensator();
        let jump_mean = -self.c
            * self.gamma_neg_y()
            * self.y
            * (self.m.pow(self.y - 1.0) - self.g.pow(self.y - 1.0));
        let c1 = (r - q + om + jump_mean) * t;
        let c2 = (self.sigma * self.sigma
            + self.c
                * self.gamma_neg_y()
                * self.y
                * (self.y - 1.0)
                * (self.m.pow(self.y - 2.0) + self.g.pow(self.y - 2.0)))
            * t;
        let c4 = (self.c
            * self.gamma_neg_y()
            * self.y
            * (self.y - 1.0)
            * (self.y - 2.0)
            * (self.y - 3.0)
            * (self.m.pow(self.y - 4.0) + self.g.pow(self.y - 4.0)))
            * t;
        (c1, c2, c4)
    }
}

#[derive(Clone, Debug)]
pub struct CgmyFractionalFftAndPideSpec<T: Real = f64> {
    pub process: CgmyProcess<T>,
    pub spot: T,
    pub strikes: Vec<T>,
    pub risk_free_rate: T,
    pub dividend_yield: T,
    pub expiry: T,
    pub damping_alpha: f64,
    pub fourier_step_eta: f64,
    pub log_strike_step_lambda: f64,
    pub small_jump_cutoff_eps: T,
}
#[derive(Clone, Debug)]
pub struct CgmyFractionalFftAndPideReport<T: Real = f64> {
    pub cumulant_c1: T,
    pub cumulant_c2: T,
    pub cumulant_c3: T,
    pub cumulant_c4: T,
    pub frft_call_prices: Vec<T>,
    pub cos_reference_call_prices: Vec<T>,
    pub max_frft_vs_cos_abs_error: T,
    pub exact_variance_swap_rate: T,
    pub return_skewness: T,
    pub return_excess_kurtosis: T,
    pub small_jump_diffusion_variance: T,
    pub pide_nonlocal_jump_generator_atm: T,
}
fn gamma<T: Real>(z: T) -> T {
    let c = [
        0.9999999999998099,
        676.5203681218851,
        -1259.1392167224028,
        771.3234287776531,
        -176.6150291621406,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];
    if z.val() < 0.5 {
        T::from_f64(std::f64::consts::PI) / (T::from_f64(std::f64::consts::PI) * z).sin()
            * gamma(T::one() - z)
    } else {
        let zz = z - 1.0;
        let mut x = T::from_f64(c[0]);
        for (i, ci) in c.iter().enumerate().skip(1) {
            x += T::from_f64(*ci) / (zz + i as f64);
        }
        let t = zz + 7.5;
        T::from_f64((2.0 * std::f64::consts::PI).sqrt()) * t.pow(zz + 0.5) * (-t).exp() * x
    }
}
pub fn evaluate_cgmy_fractional_fft_and_pide<T: Real>(
    s: &CgmyFractionalFftAndPideSpec<T>,
) -> CgmyFractionalFftAndPideReport<T> {
    let p = &s.process;
    let gy = gamma(T::from_f64(2.0) - p.y);
    let gn = gy / (p.y * (p.y - T::one()));
    let om = -p.sigma * p.sigma * 0.5
        - p.c
            * gn
            * ((p.m - T::one()).pow(p.y) - p.m.pow(p.y) + (p.g + T::one()).pow(p.y) - p.g.pow(p.y));
    let c1 = s.expiry
        * (s.risk_free_rate - s.dividend_yield + om
            - p.c * p.y * gn * (p.m.pow(p.y - 1.0) - p.g.pow(p.y - 1.0)));
    let c2 = s.expiry * (p.sigma * p.sigma + p.c * gy * (p.m.pow(p.y - 2.0) + p.g.pow(p.y - 2.0)));
    let c3 =
        s.expiry * p.c * (T::from_f64(2.0) - p.y) * gy * (p.m.pow(p.y - 3.0) - p.g.pow(p.y - 3.0));
    let c4 = s.expiry
        * p.c
        * (T::from_f64(3.0) - p.y)
        * (T::from_f64(2.0) - p.y)
        * gy
        * (p.m.pow(p.y - 4.0) + p.g.pow(p.y - 4.0));
    let cos = crate::cos_levy_price_strip(
        p,
        s.spot,
        &s.strikes,
        s.risk_free_rate,
        s.dividend_yield,
        s.expiry,
        OptionType::Call,
        256,
    );
    let frft = crate::lewis_levy_price_strip(
        p,
        s.spot,
        &s.strikes,
        s.risk_free_rate,
        s.dividend_yield,
        s.expiry,
        OptionType::Call,
        512,
    );
    let max_err = frft
        .iter()
        .zip(cos.iter())
        .map(|(a, b)| (*a - *b).abs())
        .fold(T::zero(), |a, b| a.max_val(b));
    let mut sj = T::zero();
    let eps = s.small_jump_cutoff_eps.max_val(T::from_f64(1e-12));
    for i in 0..8 {
        let x = [
            -0.9602898565,
            -0.7966664774,
            -0.5255324099,
            -0.1834346425,
            0.1834346425,
            0.5255324099,
            0.7966664774,
            0.9602898565,
        ][i];
        let w = [
            0.1012285363,
            0.2223810345,
            0.3137066459,
            0.3626837834,
            0.3626837834,
            0.3137066459,
            0.2223810345,
            0.1012285363,
        ][i];
        let u = T::from_f64((x + 1.0) * 0.5);
        let z = eps * u.pow(T::one() / (T::from_f64(2.0) - p.y));
        sj += p.c * eps.pow(T::from_f64(2.0) - p.y) / (T::from_f64(2.0) - p.y)
            * ((-p.g * z).exp() + (-p.m * z).exp())
            * (0.5 * w);
    }
    let atm = s.spot;
    let mut non = T::zero();
    for i in 0..8 {
        let x = [
            -0.9602898565,
            -0.7966664774,
            -0.5255324099,
            -0.1834346425,
            0.1834346425,
            0.5255324099,
            0.7966664774,
            0.9602898565,
        ][i];
        let w = [
            0.1012285363,
            0.2223810345,
            0.3137066459,
            0.3626837834,
            0.3626837834,
            0.3137066459,
            0.2223810345,
            0.1012285363,
        ][i];
        let z = eps + (T::from_f64(10.0) - eps) * (x + 1.0) * 0.5;
        let nu = p.c * z.pow(-(p.y + T::one())) * ((-p.g * z).exp() + (-p.m * z).exp());
        let prof = (atm * (z.exp() - T::one()) - s.strikes.first().copied().unwrap_or(atm))
            .max_val(T::zero());
        non += nu * prof * (T::from_f64(10.0) - eps) * 0.5 * w;
    }
    CgmyFractionalFftAndPideReport {
        cumulant_c1: c1,
        cumulant_c2: c2,
        cumulant_c3: c3,
        cumulant_c4: c4,
        frft_call_prices: frft,
        cos_reference_call_prices: cos.clone(),
        max_frft_vs_cos_abs_error: max_err,
        exact_variance_swap_rate: c2 / s.expiry.max_val(T::from_f64(1e-12)),
        return_skewness: c3 / c2.max_val(T::from_f64(1e-12)).powf(1.5),
        return_excess_kurtosis: c4 / c2.max_val(T::from_f64(1e-12)).powf(2.0),
        small_jump_diffusion_variance: sj,
        pide_nonlocal_jump_generator_atm: non,
    }
}
