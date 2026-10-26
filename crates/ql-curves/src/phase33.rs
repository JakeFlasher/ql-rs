use ql_ad::Real;
#[derive(Clone, Debug)]
pub struct BrigoMercurioRebonatoCmsSpec<T: Real = f64> {
    pub fixing_tenors: Vec<T>,
    pub accrual_taus: Vec<T>,
    pub forward_rates: Vec<T>,
    pub discount_factors: Vec<T>,
    pub forward_vols: Vec<T>,
    pub corr_decay_beta: T,
    pub option_expiry: T,
    pub swaption_strike: T,
}
#[derive(Clone, Copy, Debug)]
pub struct BrigoMercurioRebonatoCmsReport<T: Real = f64> {
    pub swap_annuity_pv01: T,
    pub forward_swap_rate: T,
    pub rebonato_lmm_swaption_vol: T,
    pub hagan_annuity_convexity_ratio: T,
    pub cms_convexity_adjustment_bps: T,
    pub convexity_adjusted_cms_rate: T,
    pub black76_payer_swaption_pv: T,
    pub weights_sum_conservation_error: T,
}
pub fn evaluate_brigo_mercurio_rebonato_swaption_and_cms<T: Real>(
    s: &BrigoMercurioRebonatoCmsSpec<T>,
) -> BrigoMercurioRebonatoCmsReport<T> {
    let n = s.accrual_taus.len();
    let mut a = T::zero();
    for i in 0..n {
        a += s.accrual_taus[i] * s.discount_factors[i];
    }
    let mut sw = T::zero();
    for i in 0..n {
        sw += s.accrual_taus[i] * s.discount_factors[i] * s.forward_rates[i];
    }
    let sw = sw / a;
    let mut var = T::zero();
    let mut ws = T::zero();
    for i in 0..n {
        let wi = s.accrual_taus[i] * s.discount_factors[i] / a;
        ws += wi;
        for j in 0..n {
            let wj = s.accrual_taus[j] * s.discount_factors[j] / a;
            let d = if i == j {
                T::zero()
            } else {
                let z = s.fixing_tenors[i] - s.fixing_tenors[j];
                (z * z + T::from_f64(1e-30)).sqrt()
            };
            let rho = (-(s.corr_decay_beta * d)).exp();
            var += wi
                * wj
                * s.forward_rates[i]
                * s.forward_rates[j]
                * rho
                * s.forward_vols[i]
                * s.forward_vols[j]
                / (sw * sw);
        }
    }
    let vol = var.sqrt();
    let x = sw.max_val(T::from_f64(1e-8));
    let one = T::one();
    let nf = T::from_f64(n as f64);
    let u = one + x;
    let g = (one - u.powf(-(n as f64))) / x;
    let gp = (nf * u.powf(-(n as f64) - 1.0)) / x - g / x;
    let gpp = -(nf * (nf + one) * u.powf(-(n as f64) - 2.0)) / x - T::from_f64(2.0) * gp / x;
    let ratio = -(gpp / gp);
    let adj = T::from_f64(0.5) * x * x * ratio * ((vol * vol * s.option_expiry).exp() - one);
    let cms = x + adj;
    let d1 = ((sw / s.swaption_strike).ln() + T::from_f64(0.5) * vol * vol * s.option_expiry)
        / (vol * s.option_expiry.sqrt());
    let d2 = d1 - vol * s.option_expiry.sqrt();
    let pv = a * (sw * d1.norm_cdf() - s.swaption_strike * d2.norm_cdf());
    BrigoMercurioRebonatoCmsReport {
        swap_annuity_pv01: a,
        forward_swap_rate: sw,
        rebonato_lmm_swaption_vol: vol,
        hagan_annuity_convexity_ratio: ratio,
        cms_convexity_adjustment_bps: adj * T::from_f64(10000.0),
        convexity_adjusted_cms_rate: cms,
        black76_payer_swaption_pv: pv,
        weights_sum_conservation_error: (ws - one).abs(),
    }
}
#[cfg(test)]
mod phase33_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> BrigoMercurioRebonatoCmsSpec<T> {
        BrigoMercurioRebonatoCmsSpec {
            fixing_tenors: vec![T::one(), T::from_f64(2.)],
            accrual_taus: vec![T::one(), T::one()],
            forward_rates: vec![T::from_f64(0.03) + x, T::from_f64(0.035)],
            discount_factors: vec![T::from_f64(0.98), T::from_f64(0.95)],
            forward_vols: vec![T::from_f64(0.2), T::from_f64(0.2)],
            corr_decay_beta: T::from_f64(0.1),
            option_expiry: T::one(),
            swaption_strike: T::from_f64(0.03),
        }
    }
    #[test]
    fn phase33_rebonato_smoke() {
        let r = evaluate_brigo_mercurio_rebonato_swaption_and_cms(&spec(0.));
        assert!(r.weights_sum_conservation_error < 1e-12 && r.cms_convexity_adjustment_bps > 0.);
        let h = 1e-6;
        let fd = (evaluate_brigo_mercurio_rebonato_swaption_and_cms(&spec(h))
            .black76_payer_swaption_pv
            - evaluate_brigo_mercurio_rebonato_swaption_and_cms(&spec(-h))
                .black76_payer_swaption_pv)
            / (2. * h);
        let ad =
            evaluate_brigo_mercurio_rebonato_swaption_and_cms(&spec(Dual::<1>::variable(0., 0)))
                .black76_payer_swaption_pv;
        assert!((ad.d[0] - fd).abs() < 1e-5);
    }
}
