use ql_ad::{Dual, Real, TapeScope, TapeVar};

#[derive(Clone, Copy, Debug)]
pub struct SeasonalParams<T: Real> {
    pub a1: T,
    pub b1: T,
    pub a2: T,
    pub b2: T,
}
impl<T: Real> SeasonalParams<T> {
    pub fn factor(&self, t: T) -> T {
        self.a1 * (T::from_f64(2.0 * std::f64::consts::PI) * t).cos()
            + self.b1 * (T::from_f64(2.0 * std::f64::consts::PI) * t).sin()
            + self.a2 * (T::from_f64(4.0 * std::f64::consts::PI) * t).cos()
            + self.b2 * (T::from_f64(4.0 * std::f64::consts::PI) * t).sin()
    }
}

#[derive(Clone, Debug)]
pub struct CommodityFuturesCurve<T: Real> {
    pub times: Vec<T>,
    pub log_base_forwards: Vec<T>,
    pub seasonal: SeasonalParams<T>,
}
impl<T: Real> CommodityFuturesCurve<T> {
    /// Construct a curve from futures pillars with a flat (zero) seasonal
    /// component.  The convenience constructor keeps storage/arbitrage
    /// callers concise while `from_futures_pillars` remains available for
    /// calibrated seasonal curves.
    pub fn new(times: Vec<f64>, market_futures: Vec<T>) -> Result<Self, String> {
        if times.is_empty()
            || times.len() != market_futures.len()
            || times.iter().any(|&t| !t.is_finite() || t < 0.0)
            || market_futures
                .iter()
                .any(|f| !f.val().is_finite() || f.val() <= 0.0)
        {
            return Err("invalid futures curve pillars".to_string());
        }
        let times_t = times.into_iter().map(T::from_f64).collect::<Vec<_>>();
        Ok(Self::from_futures_pillars(
            times_t,
            &market_futures,
            SeasonalParams {
                a1: T::zero(),
                b1: T::zero(),
                a2: T::zero(),
                b2: T::zero(),
            },
        ))
    }

    pub fn from_futures_pillars(
        times: Vec<T>,
        market_futures: &[T],
        seasonal: SeasonalParams<T>,
    ) -> Self {
        assert_eq!(times.len(), market_futures.len());
        assert!(!times.is_empty());
        let log_base_forwards = times
            .iter()
            .zip(market_futures)
            .map(|(t, f)| f.ln() - seasonal.factor(*t))
            .collect();
        Self {
            times,
            log_base_forwards,
            seasonal,
        }
    }
    pub fn futures_price(&self, t: T) -> T {
        let n = self.times.len();
        if t.val() <= self.times[0].val() {
            return (self.log_base_forwards[0] + self.seasonal.factor(t)).exp();
        }
        let i = if t.val() <= self.times[0].val() {
            0
        } else if t.val() >= self.times[n - 1].val() {
            n - 1
        } else {
            self.times
                .partition_point(|x| x.val() <= t.val())
                .saturating_sub(1)
        };
        let base = if i + 1 >= n {
            self.log_base_forwards[i]
        } else {
            let w = (t - self.times[i]) / (self.times[i + 1] - self.times[i]);
            self.log_base_forwards[i]
                + w * (self.log_base_forwards[i + 1] - self.log_base_forwards[i])
        };
        (base + self.seasonal.factor(t)).exp()
    }
    pub fn implied_convenience_yield(&self, spot: T, rate: T, storage_cost: T, t: T) -> T {
        rate + storage_cost - (self.futures_price(t) / spot).ln() / t
    }
}

#[derive(Clone, Debug)]
pub struct CalibratedCommodityCurve {
    pub times: Vec<f64>,
    pub market_futures: Vec<f64>,
    pub log_base_forwards: Vec<f64>,
    pub seasonal: SeasonalParams<f64>,
}
impl CalibratedCommodityCurve {
    pub fn new(times: Vec<f64>, market_futures: Vec<f64>, seasonal: SeasonalParams<f64>) -> Self {
        assert_eq!(times.len(), market_futures.len());
        let log_base_forwards = times
            .iter()
            .zip(&market_futures)
            .map(|(&t, &f)| f.ln() - seasonal.factor(t))
            .collect();
        Self {
            times,
            market_futures,
            log_base_forwards,
            seasonal,
        }
    }
    pub fn curve(&self) -> CommodityFuturesCurve<f64> {
        CommodityFuturesCurve::from_futures_pillars(
            self.times.clone(),
            &self.market_futures,
            self.seasonal,
        )
    }
    pub fn dp_dfutures<F>(&self, pricer: F) -> (f64, Vec<f64>)
    where
        F: Fn(&CommodityFuturesCurve<Dual<1>>) -> Dual<1>,
    {
        let mut out = Vec::with_capacity(self.market_futures.len());
        for j in 0..self.market_futures.len() {
            let times = self
                .times
                .iter()
                .map(|&x| Dual::from_value(x))
                .collect::<Vec<_>>();
            let theta = self
                .log_base_forwards
                .iter()
                .enumerate()
                .map(|(i, &x)| {
                    if i == j {
                        Dual::variable(x, 0)
                    } else {
                        Dual::from_value(x)
                    }
                })
                .collect::<Vec<_>>();
            let seasonal = SeasonalParams {
                a1: Dual::from_value(self.seasonal.a1),
                b1: Dual::from_value(self.seasonal.b1),
                a2: Dual::from_value(self.seasonal.a2),
                b2: Dual::from_value(self.seasonal.b2),
            };
            out.push(
                pricer(&CommodityFuturesCurve {
                    times,
                    log_base_forwards: theta,
                    seasonal,
                })
                .d[0]
                    / self.market_futures[j],
            );
        }
        let curve = CommodityFuturesCurve {
            times: self.times.iter().map(|&x| Dual::from_value(x)).collect(),
            log_base_forwards: self
                .log_base_forwards
                .iter()
                .map(|&x| Dual::from_value(x))
                .collect(),
            seasonal: SeasonalParams {
                a1: Dual::from_value(self.seasonal.a1),
                b1: Dual::from_value(self.seasonal.b1),
                a2: Dual::from_value(self.seasonal.a2),
                b2: Dual::from_value(self.seasonal.b2),
            },
        };
        (pricer(&curve).val, out)
    }

    /// Evaluate a portfolio once on the log-forward state tape and lift the
    /// adjoint through `d log(F_j) / d F_j = 1/F_j`.  The seasonal factors are
    /// part of the curve state, not of the market quote, so no second
    /// seasonal subtraction occurs in this VJP.
    pub fn value_and_par_futures_delta<F>(&self, pricer: F) -> (f64, Vec<f64>)
    where
        F: FnOnce(&CommodityFuturesCurve<TapeVar<f64>>) -> TapeVar<f64>,
    {
        let guard = TapeScope::guard();
        let vars: Vec<_> = self
            .log_base_forwards
            .iter()
            .map(|&x| TapeScope::with_active(|t| t.input(x)))
            .collect();
        let curve = CommodityFuturesCurve {
            times: self.times.iter().map(|&x| TapeVar::from_f64(x)).collect(),
            log_base_forwards: vars.clone(),
            seasonal: SeasonalParams {
                a1: TapeVar::from_f64(self.seasonal.a1),
                b1: TapeVar::from_f64(self.seasonal.b1),
                a2: TapeVar::from_f64(self.seasonal.a2),
                b2: TapeVar::from_f64(self.seasonal.b2),
            },
        };
        let out = pricer(&curve);
        let adj = TapeScope::gradient_from(out);
        let deltas = vars
            .iter()
            .enumerate()
            .map(|(j, v)| adj.gradient(*v) / self.market_futures[j])
            .collect();
        let value = out.val;
        let _ = guard.stop();
        (value, deltas)
    }
}

mod evaluate_routledge_seppi_spatt_stockout_backwardation_module {
    // Routledge, Seppi & Spatt (2000, JF) & Geman (2005) Competitive Storage Rational
    // Expectations Equilibrium, Stockout Boundary & Endogenous Backwardation Volatility Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct RoutledgeSeppiSpattStockoutBackwardationSpec<T: Real = f64> {
        pub spot_commodity_price_s0: T,
        pub discretionary_inventory_i: T,
        pub critical_stockout_threshold_i_star: T,
        pub base_abundant_volatility_sigma0: T,
        pub stockout_volatility_spike_eta: T,
        pub interest_plus_storage_cost_rate_c: T,
        pub convenience_scarcity_slope_alpha: T,
        pub maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct RoutledgeSeppiSpattStockoutBackwardationReport<T: Real = f64> {
        pub inventory_cushion_above_stockout: T,
        pub endogenous_spot_volatility: T,
        pub endogenous_convenience_yield: T,
        pub spot_forward_correlation_rho: T,
        pub equilibrium_forward_price: T,
        pub backwardation_basis_pct: T,
    }

    pub fn evaluate_routledge_seppi_spatt_stockout_backwardation<T: Real>(
        spec: &RoutledgeSeppiSpattStockoutBackwardationSpec<T>,
    ) -> RoutledgeSeppiSpattStockoutBackwardationReport<T> {
        let s0 = spec.spot_commodity_price_s0.max_val(T::from_f64(1.0));
        let inv = spec.discretionary_inventory_i.max_val(T::from_f64(0.05));
        let i_star = spec
            .critical_stockout_threshold_i_star
            .max_val(T::from_f64(0.01));
        let sig0 = spec
            .base_abundant_volatility_sigma0
            .max_val(T::from_f64(0.02));
        let eta = spec.stockout_volatility_spike_eta.max_val(T::zero());
        let cost_c = spec.interest_plus_storage_cost_rate_c;
        let alpha = spec.convenience_scarcity_slope_alpha.max_val(T::zero());
        let t = spec.maturity_years_t.max_val(T::from_f64(0.05));

        let raw_gap = inv - i_star;
        let eps = T::from_f64(0.05);
        let cushion = T::from_f64(0.5) * (raw_gap + (raw_gap * raw_gap + eps * eps).sqrt());

        // Routledge-Seppi-Spatt (2000): volatility and convenience yield surge as inventory approaches I*:
        let scarcity_kernel = (-cushion / i_star.max_val(T::from_f64(0.10))).exp();
        let sig_spot = sig0 * (T::one() + eta * scarcity_kernel);

        // Endogenous spot-forward correlation drops below 1 during stockouts because short-term spikes mean-revert:
        let rho_sf = T::one() / (T::one() + T::from_f64(0.35) * eta * scarcity_kernel * t).sqrt();

        // Routledge, Seppi & Spatt (2000) endogenous convenience yield includes timing-option variance premium:
        let conv_y = (alpha
            + T::from_f64(0.5) * sig_spot * sig_spot * (T::one() - rho_sf * rho_sf))
            * scarcity_kernel;
        let fwd = s0 * ((cost_c - conv_y) * t).exp();
        let back_pct = ((s0 - fwd) / s0) * T::from_f64(100.0);

        RoutledgeSeppiSpattStockoutBackwardationReport {
            inventory_cushion_above_stockout: cushion,
            endogenous_spot_volatility: sig_spot,
            endogenous_convenience_yield: conv_y,
            spot_forward_correlation_rho: rho_sf,
            equilibrium_forward_price: fwd,
            backwardation_basis_pct: back_pct,
        }
    }

    #[cfg(test)]
    mod evaluate_routledge_seppi_spatt_stockout_backwardation_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn routledge_seppi_spatt_stockout_backwardation_invariants_and_ad_parity() {
            let spec = RoutledgeSeppiSpattStockoutBackwardationSpec {
                spot_commodity_price_s0: 78.0,
                discretionary_inventory_i: 0.24,
                critical_stockout_threshold_i_star: 0.20,
                base_abundant_volatility_sigma0: 0.22,
                stockout_volatility_spike_eta: 0.85,
                interest_plus_storage_cost_rate_c: 0.065,
                convenience_scarcity_slope_alpha: 0.16,
                maturity_years_t: 0.50,
            };
            let rep = evaluate_routledge_seppi_spatt_stockout_backwardation(&spec);
            assert!(rep.endogenous_spot_volatility > spec.base_abundant_volatility_sigma0);
            assert!(rep.backwardation_basis_pct > 0.0);
            assert!(rep.spot_forward_correlation_rho < 1.0);

            let h = 1e-6;
            let r_ad = evaluate_routledge_seppi_spatt_stockout_backwardation(
                &RoutledgeSeppiSpattStockoutBackwardationSpec {
                    spot_commodity_price_s0: Dual::<8>::variable(spec.spot_commodity_price_s0, 0),
                    discretionary_inventory_i: Dual::<8>::variable(
                        spec.discretionary_inventory_i,
                        1,
                    ),
                    critical_stockout_threshold_i_star: Dual::<8>::variable(
                        spec.critical_stockout_threshold_i_star,
                        2,
                    ),
                    base_abundant_volatility_sigma0: Dual::<8>::variable(
                        spec.base_abundant_volatility_sigma0,
                        3,
                    ),
                    stockout_volatility_spike_eta: Dual::<8>::variable(
                        spec.stockout_volatility_spike_eta,
                        4,
                    ),
                    interest_plus_storage_cost_rate_c: Dual::<8>::variable(
                        spec.interest_plus_storage_cost_rate_c,
                        5,
                    ),
                    convenience_scarcity_slope_alpha: Dual::<8>::variable(
                        spec.convenience_scarcity_slope_alpha,
                        6,
                    ),
                    maturity_years_t: Dual::<8>::variable(spec.maturity_years_t, 7),
                },
            );
            for (idx, &ad) in r_ad.equilibrium_forward_price.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_commodity_price_s0 += h;
                        sm.spot_commodity_price_s0 -= h;
                    }
                    1 => {
                        sp.discretionary_inventory_i += h;
                        sm.discretionary_inventory_i -= h;
                    }
                    2 => {
                        sp.critical_stockout_threshold_i_star += h;
                        sm.critical_stockout_threshold_i_star -= h;
                    }
                    3 => {
                        sp.base_abundant_volatility_sigma0 += h;
                        sm.base_abundant_volatility_sigma0 -= h;
                    }
                    4 => {
                        sp.stockout_volatility_spike_eta += h;
                        sm.stockout_volatility_spike_eta -= h;
                    }
                    5 => {
                        sp.interest_plus_storage_cost_rate_c += h;
                        sm.interest_plus_storage_cost_rate_c -= h;
                    }
                    6 => {
                        sp.convenience_scarcity_slope_alpha += h;
                        sm.convenience_scarcity_slope_alpha -= h;
                    }
                    _ => {
                        sp.maturity_years_t += h;
                        sm.maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_routledge_seppi_spatt_stockout_backwardation(&sp)
                    .equilibrium_forward_price
                    - evaluate_routledge_seppi_spatt_stockout_backwardation(&sm)
                        .equilibrium_forward_price)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_routledge_seppi_spatt_stockout_backwardation_module::*;
