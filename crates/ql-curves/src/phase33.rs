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
    let mut wsum = T::zero();
    for i in 0..n {
        let wi = s.accrual_taus[i] * s.discount_factors[i] / a;
        wsum += wi;
        for j in 0..n {
            let wj = s.accrual_taus[j] * s.discount_factors[j] / a;
            let dt = (s.fixing_tenors[i] - s.fixing_tenors[j]).abs();
            let rho = (-(s.corr_decay_beta * dt)).exp();
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
    let x = sw;
    let one = T::one();
    let g = (one - (one + x).powf(-(n as f64))) / x;
    let h = T::from_f64(1e-5);
    let gp = ((one - (one + x + h).powf(-(n as f64))) / (x + h)
        - (one - (one + x - h).powf(-(n as f64))) / (x - h))
        / (h * 2.0);
    let gpp = (((one - (one + x + h).powf(-(n as f64))) / (x + h)) - T::from_f64(2.0) * g
        + ((one - (one + x - h).powf(-(n as f64))) / (x - h)))
        / (h * h);
    let ratio = -(gpp / gp);
    let adj = T::from_f64(0.5) * x * x * ratio * ((vol * vol * s.option_expiry).exp() - one);
    let cms = x + adj;
    let d1 = ((cms / s.swaption_strike).ln() + T::from_f64(0.5) * vol * vol * s.option_expiry)
        / (vol * s.option_expiry.sqrt());
    let d2 = d1 - vol * s.option_expiry.sqrt();
    let pv = a * (cms * d1.norm_cdf() - s.swaption_strike * d2.norm_cdf());
    BrigoMercurioRebonatoCmsReport {
        swap_annuity_pv01: a,
        forward_swap_rate: sw,
        rebonato_lmm_swaption_vol: vol,
        hagan_annuity_convexity_ratio: ratio,
        cms_convexity_adjustment_bps: adj * T::from_f64(10000.0),
        convexity_adjusted_cms_rate: cms,
        black76_payer_swaption_pv: pv,
        weights_sum_conservation_error: (wsum - one).abs(),
    }
}
#[cfg(test)]
mod phase33_tests {
    use super::*;
    #[test]
    fn phase33_rebonato_smoke() {
        let s = BrigoMercurioRebonatoCmsSpec {
            fixing_tenors: vec![1.0, 2.0],
            accrual_taus: vec![1.0, 1.0],
            forward_rates: vec![0.03, 0.035],
            discount_factors: vec![0.98, 0.95],
            forward_vols: vec![0.2, 0.2],
            corr_decay_beta: 0.1,
            option_expiry: 1.0,
            swaption_strike: 0.03,
        };
        assert!(
            evaluate_brigo_mercurio_rebonato_swaption_and_cms(&s).rebonato_lmm_swaption_vol > 0.0
        );
    }
}
