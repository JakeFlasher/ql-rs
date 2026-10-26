use ql_ad::Real;
use ql_math::{inv_norm_cdf, norm_cdf};

#[derive(Clone, Copy, Debug)]
pub struct ClarkWystupFxVannaVolgaSmileSpec<T: Real = f64> {
    pub spot_fx: T,
    pub domestic_rate: T,
    pub foreign_rate: T,
    pub expiry_years: T,
    pub atm_vol: T,
    pub risk_reversal_25d: T,
    pub butterfly_25d: T,
    pub target_strike: T,
}

#[derive(Clone, Copy, Debug)]
pub struct ClarkWystupFxVannaVolgaSmileReport<T: Real = f64> {
    pub vol_25d_put: T,
    pub vol_25d_call: T,
    pub pillar_reconstruction_error: T,
    pub strike_25d_put: T,
    pub strike_atm_dns: T,
    pub strike_25d_call: T,
    pub vanna_volga_weights: [T; 3],
    pub weights_partition_of_unity_error: T,
    pub vanna_volga_implied_vol: T,
    pub vanna_volga_call_price: T,
    pub flat_atm_bs_call_price: T,
    pub smile_overhedge_adjustment: T,
}

fn gk_call<T: Real>(s: T, k: T, rd: T, rf: T, sigma: T, t: T) -> T {
    let st = sigma * t.sqrt();
    let d1 = ((s / k).ln() + (rd - rf + T::from_f64(0.5) * sigma * sigma) * t) / st;
    let d2 = d1 - st;
    s * (-rf * t).exp() * norm_cdf(d1) - k * (-rd * t).exp() * norm_cdf(d2)
}

pub fn evaluate_clark_wystup_fx_vanna_volga_smile<T: Real>(
    s: &ClarkWystupFxVannaVolgaSmileSpec<T>,
) -> ClarkWystupFxVannaVolgaSmileReport<T> {
    let t = s.expiry_years.max_val(T::from_f64(1e-6));
    let half = T::from_f64(0.5);
    let sig_c = s.atm_vol + s.butterfly_25d + half * s.risk_reversal_25d;
    let sig_p = s.atm_vol + s.butterfly_25d - half * s.risk_reversal_25d;
    let rec1 = (sig_c - sig_p) - s.risk_reversal_25d;
    let rec2 = half * (sig_c + sig_p) - s.atm_vol - s.butterfly_25d;
    let rec_err = (rec1 * rec1 + rec2 * rec2 + T::from_f64(1e-30)).sqrt();
    let fwd = s.spot_fx * ((s.domestic_rate - s.foreign_rate) * t).exp();
    let alpha_d = inv_norm_cdf(
        (T::from_f64(0.25) * (s.foreign_rate * t).exp())
            .max_val(T::from_f64(0.01))
            .min_val(T::from_f64(0.49)),
    );
    let k_atm = fwd * (half * s.atm_vol * s.atm_vol * t).exp();
    let k1 = fwd * (alpha_d * sig_p * t.sqrt() + half * sig_p * sig_p * t).exp();
    let k3 = fwd * (-alpha_d * sig_c * t.sqrt() + half * sig_c * sig_c * t).exp();
    let x1 = k1.ln();
    let x2 = k_atm.ln();
    let x3 = k3.ln();
    let x = s.target_strike.ln();
    let w1 = (x - x2) * (x - x3) / ((x1 - x2) * (x1 - x3));
    let w2 = (x - x1) * (x - x3) / ((x2 - x1) * (x2 - x3));
    let w3 = (x - x1) * (x - x2) / ((x3 - x1) * (x3 - x2));
    let w_diff = T::one() - (w1 + w2 + w3);
    let w_err = (w_diff * w_diff + T::from_f64(1e-30)).sqrt();
    let sig_vv = w1 * sig_p + w2 * s.atm_vol + w3 * sig_c;
    let flat_c = gk_call(
        s.spot_fx,
        s.target_strike,
        s.domestic_rate,
        s.foreign_rate,
        s.atm_vol,
        t,
    );
    let vv_c = gk_call(
        s.spot_fx,
        s.target_strike,
        s.domestic_rate,
        s.foreign_rate,
        sig_vv,
        t,
    );
    ClarkWystupFxVannaVolgaSmileReport {
        vol_25d_put: sig_p,
        vol_25d_call: sig_c,
        pillar_reconstruction_error: rec_err,
        strike_25d_put: k1,
        strike_atm_dns: k_atm,
        strike_25d_call: k3,
        vanna_volga_weights: [w1, w2, w3],
        weights_partition_of_unity_error: w_err,
        vanna_volga_implied_vol: sig_vv,
        vanna_volga_call_price: vv_c,
        flat_atm_bs_call_price: flat_c,
        smile_overhedge_adjustment: vv_c - flat_c,
    }
}

#[cfg(test)]
mod phase35_tests {
    use super::*;
    use ql_ad::Dual;
    fn spec<T: Real>(x: T) -> ClarkWystupFxVannaVolgaSmileSpec<T> {
        ClarkWystupFxVannaVolgaSmileSpec {
            spot_fx: T::from_f64(1.10) + x,
            domestic_rate: T::from_f64(0.03),
            foreign_rate: T::from_f64(0.015),
            expiry_years: T::from_f64(1.0),
            atm_vol: T::from_f64(0.10),
            risk_reversal_25d: T::from_f64(0.012),
            butterfly_25d: T::from_f64(0.004),
            target_strike: T::from_f64(1.14),
        }
    }
    #[test]
    fn phase35_fx_invariants_and_dual_fd() {
        let r = evaluate_clark_wystup_fx_vanna_volga_smile(&spec(0.0));
        assert!(r.pillar_reconstruction_error < 1e-12);
        assert!(r.weights_partition_of_unity_error < 1e-12);
        assert!(r.strike_25d_put < r.strike_atm_dns && r.strike_atm_dns < r.strike_25d_call);
        let h = 1e-6;
        let fd = (evaluate_clark_wystup_fx_vanna_volga_smile(&spec(h)).vanna_volga_call_price
            - evaluate_clark_wystup_fx_vanna_volga_smile(&spec(-h)).vanna_volga_call_price)
            / (2.0 * h);
        let ad = evaluate_clark_wystup_fx_vanna_volga_smile(&spec(Dual::<1>::variable(0.0, 0)))
            .vanna_volga_call_price;
        assert!(ad.d[0].abs() > 1e-6 && (ad.d[0] - fd).abs() < 1e-5);
    }
}
