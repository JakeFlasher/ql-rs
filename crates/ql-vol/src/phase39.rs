//! Gatheral (2006) & Bergomi (2016) VIX Futures Concavity Discount,
//! Variance-Swap vs. VIX Payoff Gap & VVIX Log-Contract Smile Engine.

use ql_ad::Real;
use ql_math::norm_cdf;

#[derive(Debug, Clone, Copy)]
pub struct GatheralBergomiVixVvixSpec<T: Real = f64> {
    pub forward_variance_swap_var: T,
    pub vvix_vol_of_vol_omega: T,
    pub vix_option_strike_k: T,
    pub risk_free_rate: T,
    pub expiry_years_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct GatheralBergomiVixVvixReport<T: Real = f64> {
    pub fair_volatility_sqrt_var: T,
    pub vix_futures_price: T,
    pub jensen_concavity_discount: T,
    pub vix_variance_swap_replication_variance: T,
    pub black76_vix_call_price: T,
    pub vvix_log_contract_implied_vol: T,
}

pub fn evaluate_gatheral_bergomi_vix_futures_and_vvix_smile<T: Real>(
    spec: &GatheralBergomiVixVvixSpec<T>,
) -> GatheralBergomiVixVvixReport<T> {
    let v0 = spec.forward_variance_swap_var.max_val(T::from_f64(1e-6));
    let t = spec.expiry_years_t.max_val(T::from_f64(1e-3));
    let omega = spec.vvix_vol_of_vol_omega.max_val(T::from_f64(1e-4));

    // Fair forward volatility from forward variance swap: sqrt(E_0[VIX_T^2])
    let sqrt_v0 = v0.sqrt();
    // Under lognormal VIX_T dynamics with vol-of-vol omega (VVIX), E_0[VIX_T^2] = F_VIX^2 * exp(omega^2 * T)
    // Therefore F_VIX(0,T) = E_0[VIX_T] = sqrt(V_0) * exp(-0.5 * omega^2 * T)
    let f_vix = sqrt_v0 * (-T::from_f64(0.5) * omega * omega * t).exp();
    let concavity_discount = sqrt_v0 - f_vix;
    let var_vix = v0 - f_vix * f_vix;

    // Black-76 VIX call option on VIX futures F_VIX with strike K and vol-of-vol omega:
    let k = spec.vix_option_strike_k.max_val(T::from_f64(1e-4));
    let std = omega * t.sqrt();
    let d1 = ((f_vix / k).ln() + T::from_f64(0.5) * std * std) / std;
    let d2 = d1 - std;
    let df = (-spec.risk_free_rate * t).exp();
    let vix_call = df * (f_vix * norm_cdf(d1) - k * norm_cdf(d2));

    // VVIX log-contract implied volatility identity: sqrt(2/T * ln(sqrt(V_0) / F_VIX)) == omega
    let vvix_implied = ((T::from_f64(2.0) / t) * (sqrt_v0 / f_vix).ln())
        .max_val(T::zero())
        .sqrt();

    GatheralBergomiVixVvixReport {
        fair_volatility_sqrt_var: sqrt_v0,
        vix_futures_price: f_vix,
        jensen_concavity_discount: concavity_discount,
        vix_variance_swap_replication_variance: var_vix,
        black76_vix_call_price: vix_call,
        vvix_log_contract_implied_vol: vvix_implied,
    }
}

#[cfg(test)]
mod phase39_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn gatheral_bergomi_vix_vvix_invariants_and_ad_parity() {
        let spec = GatheralBergomiVixVvixSpec {
            forward_variance_swap_var: 0.0484, // 22.0 vol points squared
            vvix_vol_of_vol_omega: 0.85,
            vix_option_strike_k: 0.21,
            risk_free_rate: 0.04,
            expiry_years_t: 0.25,
        };
        let rep = evaluate_gatheral_bergomi_vix_futures_and_vvix_smile(&spec);
        // Jensen concavity: E[VIX_T] < sqrt(E[VIX_T^2])
        assert!(rep.jensen_concavity_discount > 0.0);
        assert!(rep.vix_variance_swap_replication_variance > 0.0);
        assert!((rep.vvix_log_contract_implied_vol - spec.vvix_vol_of_vol_omega).abs() < 1e-12);

        let h = 1e-6;
        let r_ad =
            evaluate_gatheral_bergomi_vix_futures_and_vvix_smile(&GatheralBergomiVixVvixSpec {
                forward_variance_swap_var: Dual::<1>::variable(spec.forward_variance_swap_var, 0),
                vvix_vol_of_vol_omega: Dual::<1>::from_value(spec.vvix_vol_of_vol_omega),
                vix_option_strike_k: Dual::<1>::from_value(spec.vix_option_strike_k),
                risk_free_rate: Dual::<1>::from_value(spec.risk_free_rate),
                expiry_years_t: Dual::<1>::from_value(spec.expiry_years_t),
            });
        let mut sp = spec;
        sp.forward_variance_swap_var += h;
        let mut sm = spec;
        sm.forward_variance_swap_var -= h;
        let fd = (evaluate_gatheral_bergomi_vix_futures_and_vvix_smile(&sp).black76_vix_call_price
            - evaluate_gatheral_bergomi_vix_futures_and_vvix_smile(&sm).black76_vix_call_price)
            / (2.0 * h);
        assert!((r_ad.black76_vix_call_price.d[0] - fd).abs() < 1e-6);
    }
}
