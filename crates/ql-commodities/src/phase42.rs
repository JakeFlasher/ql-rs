//! Geman (2005) & Working (1949) / Brennan (1958) Theory of Storage
//! Inventory-Scarcity Convenience Yield Curve & Stock-Out Volatility Engine.

use ql_ad::Real;

#[derive(Debug, Clone, Copy)]
pub struct GemanWorkingInventoryConvenienceSpec<T: Real = f64> {
    pub spot_commodity_price_s0: T,
    pub normalized_inventory_i_over_icrit: T,
    pub stockout_convenience_yield_y0: T,
    pub abundance_convenience_floor_yinf: T,
    pub working_curve_steepness_kappa: T,
    pub physical_storage_cost_rate_c: T,
    pub risk_free_rate_r: T,
    pub baseline_spot_vol_sigma0: T,
    pub scarcity_vol_amplification_eta: T,
    pub futures_maturity_t: T,
}

#[derive(Debug, Clone, Copy)]
pub struct GemanWorkingInventoryConvenienceReport<T: Real = f64> {
    pub marginal_convenience_yield_y: T,
    pub net_cost_of_carry_rate: T,
    pub storage_equilibrium_futures_price: T,
    pub backwardation_yield_spread: T,
    pub critical_switch_inventory_ratio: T,
    pub inventory_conditioned_spot_volatility: T,
}

pub fn evaluate_geman_working_curve_inventory_convenience_yield<T: Real>(
    spec: &GemanWorkingInventoryConvenienceSpec<T>,
) -> GemanWorkingInventoryConvenienceReport<T> {
    let s0 = spec.spot_commodity_price_s0.max_val(T::from_f64(1e-4));
    let inv = spec
        .normalized_inventory_i_over_icrit
        .max_val(T::from_f64(0.05));
    let kap = spec.working_curve_steepness_kappa.max_val(T::from_f64(0.1));
    let t = spec.futures_maturity_t.max_val(T::from_f64(0.05));

    // Working (1949) / Geman (2005, Ch. 3) convex convenience yield curve:
    // y(I) = y_inf + (y_0 - y_inf) * exp(-kappa * (I / I_crit))
    let scarcity_factor = (-kap * inv).exp();
    let y_conv = spec.abundance_convenience_floor_yinf
        + (spec.stockout_convenience_yield_y0 - spec.abundance_convenience_floor_yinf)
            * scarcity_factor;

    let gross_carry = spec.risk_free_rate_r + spec.physical_storage_cost_rate_c;
    let net_carry = gross_carry - y_conv;
    let fut_price = s0 * (net_carry * t).exp();
    let back_spread = y_conv - gross_carry;

    // Critical inventory ratio I^* / I_crit where y(I^*) = r + c (contango/backwardation boundary):
    let num = (gross_carry - spec.abundance_convenience_floor_yinf).max_val(T::from_f64(1e-4));
    let den = (spec.stockout_convenience_yield_y0 - spec.abundance_convenience_floor_yinf)
        .max_val(num + T::from_f64(1e-4));
    let i_star = -(num / den).ln() / kap;

    // Samuelson / Geman inventory-dependent spot volatility amplification under low stocks:
    let spot_vol = spec.baseline_spot_vol_sigma0
        * (T::one() + spec.scarcity_vol_amplification_eta * scarcity_factor);

    GemanWorkingInventoryConvenienceReport {
        marginal_convenience_yield_y: y_conv,
        net_cost_of_carry_rate: net_carry,
        storage_equilibrium_futures_price: fut_price,
        backwardation_yield_spread: back_spread,
        critical_switch_inventory_ratio: i_star,
        inventory_conditioned_spot_volatility: spot_vol,
    }
}

#[cfg(test)]
mod phase42_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn geman_working_curve_inventory_convenience_yield_invariants_and_ad_parity() {
        let spec = GemanWorkingInventoryConvenienceSpec {
            spot_commodity_price_s0: 82.0,
            normalized_inventory_i_over_icrit: 0.55,
            stockout_convenience_yield_y0: 0.32,
            abundance_convenience_floor_yinf: 0.01,
            working_curve_steepness_kappa: 1.8,
            physical_storage_cost_rate_c: 0.04,
            risk_free_rate_r: 0.04,
            baseline_spot_vol_sigma0: 0.28,
            scarcity_vol_amplification_eta: 0.65,
            futures_maturity_t: 1.0,
        };
        let rep = evaluate_geman_working_curve_inventory_convenience_yield(&spec);
        // Tight inventory (0.55 < I^*) produces positive backwardation (F < S_0):
        assert!(rep.backwardation_yield_spread > 0.0);
        assert!(rep.storage_equilibrium_futures_price < spec.spot_commodity_price_s0);
        assert!(rep.inventory_conditioned_spot_volatility > spec.baseline_spot_vol_sigma0);

        let h = 1e-6;
        let r_ad = evaluate_geman_working_curve_inventory_convenience_yield(
            &GemanWorkingInventoryConvenienceSpec {
                spot_commodity_price_s0: Dual::<1>::from_value(spec.spot_commodity_price_s0),
                normalized_inventory_i_over_icrit: Dual::<1>::variable(
                    spec.normalized_inventory_i_over_icrit,
                    0,
                ),
                stockout_convenience_yield_y0: Dual::<1>::from_value(
                    spec.stockout_convenience_yield_y0,
                ),
                abundance_convenience_floor_yinf: Dual::<1>::from_value(
                    spec.abundance_convenience_floor_yinf,
                ),
                working_curve_steepness_kappa: Dual::<1>::from_value(
                    spec.working_curve_steepness_kappa,
                ),
                physical_storage_cost_rate_c: Dual::<1>::from_value(
                    spec.physical_storage_cost_rate_c,
                ),
                risk_free_rate_r: Dual::<1>::from_value(spec.risk_free_rate_r),
                baseline_spot_vol_sigma0: Dual::<1>::from_value(spec.baseline_spot_vol_sigma0),
                scarcity_vol_amplification_eta: Dual::<1>::from_value(
                    spec.scarcity_vol_amplification_eta,
                ),
                futures_maturity_t: Dual::<1>::from_value(spec.futures_maturity_t),
            },
        );
        let mut sp = spec;
        sp.normalized_inventory_i_over_icrit += h;
        let mut sm = spec;
        sm.normalized_inventory_i_over_icrit -= h;
        let fd = (evaluate_geman_working_curve_inventory_convenience_yield(&sp)
            .storage_equilibrium_futures_price
            - evaluate_geman_working_curve_inventory_convenience_yield(&sm)
                .storage_equilibrium_futures_price)
            / (2.0 * h);
        assert!((r_ad.storage_equilibrium_futures_price.d[0] - fd).abs() < 1e-5);
    }
}
