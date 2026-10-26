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

mod evaluate_schwartz_one_factor_mean_reverting_log_futures_module {
    // Schwartz (1997, J. Finance Model 1) & Geman (2005, Ch. 3) Mean-Reverting Log-Price
    // Commodity Futures Curve, Samuelson Volatility Attenuation & Futures Option Engine.

    use ql_ad::Real;
    use ql_math::norm_cdf;

    #[derive(Debug, Clone, Copy)]
    pub struct SchwartzOneFactorLogFuturesSpec<T: Real = f64> {
        pub spot_commodity_price_s0: T,
        pub risk_neutral_long_run_log_price_mu_star: T,
        pub mean_reversion_speed_kappa: T,
        pub spot_log_volatility_sigma: T,
        pub option_expiry_years_te: T,
        pub futures_maturity_years_tf: T,
        pub option_strike_price_k: T,
        pub risk_free_rate_r: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct SchwartzOneFactorLogFuturesReport<T: Real = f64> {
        pub exact_log_futures_price: T,
        pub schwartz_equilibrium_futures_price: T,
        pub samuelson_terminal_futures_volatility: T,
        pub integrated_option_log_futures_variance: T,
        pub effective_black76_futures_vol: T,
        pub futures_call_option_price: T,
    }

    pub fn evaluate_schwartz_one_factor_mean_reverting_log_futures<T: Real>(
        spec: &SchwartzOneFactorLogFuturesSpec<T>,
    ) -> SchwartzOneFactorLogFuturesReport<T> {
        let s0 = spec.spot_commodity_price_s0.max_val(T::from_f64(1.0));
        let mu_star = spec.risk_neutral_long_run_log_price_mu_star;
        let kap = spec.mean_reversion_speed_kappa.max_val(T::from_f64(0.02));
        let sig = spec.spot_log_volatility_sigma.max_val(T::from_f64(0.02));
        let te = spec.option_expiry_years_te.max_val(T::from_f64(0.05));
        let tf = spec
            .futures_maturity_years_tf
            .max_val(te + T::from_f64(0.05));
        let k = spec.option_strike_price_k.max_val(T::from_f64(1.0));
        let r = spec.risk_free_rate_r;

        let x0 = s0.ln();
        let e_ktf = (-kap * tf).exp();
        let e_2ktf = (-T::from_f64(2.0) * kap * tf).exp();

        // Schwartz (1997, Eq. 8): ln F(0, T_f) = e^{-kappa T_f} ln S_0 + (1 - e^{-kappa T_f}) mu^* + (sigma^2 / (4 kappa))(1 - e^{-2 kappa T_f})
        let ln_f = e_ktf * x0
            + (T::one() - e_ktf) * mu_star
            + (sig * sig / (T::from_f64(4.0) * kap)) * (T::one() - e_2ktf);
        let f0 = ln_f.exp();

        // Samuelson effect: instantaneous volatility of F(t, T_f) at t = T_e is sigma * exp(-kappa (T_f - T_e)):
        let sam_vol = sig * (-kap * (tf - te)).exp();

        // Integrated variance of ln F(T_e, T_f) over [0, T_e] (Clewlow & Strickland 2000, Eq. 2.16):
        let var_ln_f = (sig * sig / (T::from_f64(2.0) * kap))
            * (-T::from_f64(2.0) * kap * (tf - te)).exp()
            * (T::one() - (-T::from_f64(2.0) * kap * te).exp());
        let std_f = var_ln_f.sqrt();
        let eff_vol = (var_ln_f / te).sqrt();

        let d1 = ((f0 / k).ln() + T::from_f64(0.5) * var_ln_f) / std_f;
        let d2 = d1 - std_f;
        let call_val = (-r * te).exp() * (f0 * norm_cdf(d1) - k * norm_cdf(d2));

        SchwartzOneFactorLogFuturesReport {
            exact_log_futures_price: ln_f,
            schwartz_equilibrium_futures_price: f0,
            samuelson_terminal_futures_volatility: sam_vol,
            integrated_option_log_futures_variance: var_ln_f,
            effective_black76_futures_vol: eff_vol,
            futures_call_option_price: call_val,
        }
    }

    #[cfg(test)]
    mod evaluate_schwartz_one_factor_mean_reverting_log_futures_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn schwartz_one_factor_mean_reverting_log_futures_invariants_and_ad_parity() {
            let spec = SchwartzOneFactorLogFuturesSpec {
                spot_commodity_price_s0: 76.0,
                risk_neutral_long_run_log_price_mu_star: 4.25,
                mean_reversion_speed_kappa: 0.85,
                spot_log_volatility_sigma: 0.34,
                option_expiry_years_te: 0.50,
                futures_maturity_years_tf: 1.00,
                option_strike_price_k: 74.0,
                risk_free_rate_r: 0.042,
            };
            let rep = evaluate_schwartz_one_factor_mean_reverting_log_futures(&spec);
            assert!(rep.samuelson_terminal_futures_volatility < spec.spot_log_volatility_sigma);
            assert!(rep.effective_black76_futures_vol < spec.spot_log_volatility_sigma);
            assert!(rep.futures_call_option_price > 0.0);

            let h = 1e-5;
            let r_ad = evaluate_schwartz_one_factor_mean_reverting_log_futures(
                &SchwartzOneFactorLogFuturesSpec {
                    spot_commodity_price_s0: Dual::<8>::variable(spec.spot_commodity_price_s0, 0),
                    risk_neutral_long_run_log_price_mu_star: Dual::<8>::variable(
                        spec.risk_neutral_long_run_log_price_mu_star,
                        1,
                    ),
                    mean_reversion_speed_kappa: Dual::<8>::variable(
                        spec.mean_reversion_speed_kappa,
                        2,
                    ),
                    spot_log_volatility_sigma: Dual::<8>::variable(
                        spec.spot_log_volatility_sigma,
                        3,
                    ),
                    option_expiry_years_te: Dual::<8>::variable(spec.option_expiry_years_te, 4),
                    futures_maturity_years_tf: Dual::<8>::variable(
                        spec.futures_maturity_years_tf,
                        5,
                    ),
                    option_strike_price_k: Dual::<8>::variable(spec.option_strike_price_k, 6),
                    risk_free_rate_r: Dual::<8>::variable(spec.risk_free_rate_r, 7),
                },
            );
            for (idx, &ad) in r_ad.futures_call_option_price.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_commodity_price_s0 += h;
                        sm.spot_commodity_price_s0 -= h;
                    }
                    1 => {
                        sp.risk_neutral_long_run_log_price_mu_star += h;
                        sm.risk_neutral_long_run_log_price_mu_star -= h;
                    }
                    2 => {
                        sp.mean_reversion_speed_kappa += h;
                        sm.mean_reversion_speed_kappa -= h;
                    }
                    3 => {
                        sp.spot_log_volatility_sigma += h;
                        sm.spot_log_volatility_sigma -= h;
                    }
                    4 => {
                        sp.option_expiry_years_te += h;
                        sm.option_expiry_years_te -= h;
                    }
                    5 => {
                        sp.futures_maturity_years_tf += h;
                        sm.futures_maturity_years_tf -= h;
                    }
                    6 => {
                        sp.option_strike_price_k += h;
                        sm.option_strike_price_k -= h;
                    }
                    _ => {
                        sp.risk_free_rate_r += h;
                        sm.risk_free_rate_r -= h;
                    }
                }
                let fd = (evaluate_schwartz_one_factor_mean_reverting_log_futures(&sp)
                    .futures_call_option_price
                    - evaluate_schwartz_one_factor_mean_reverting_log_futures(&sm)
                        .futures_call_option_price)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-4);
            }
        }
    }
}
pub use evaluate_schwartz_one_factor_mean_reverting_log_futures_module::*;

mod evaluate_working_brennan_storage_convenience_yield_curve_module {
    // Working (1949, AER, "The Theory of Price of Storage"), Brennan (1958, AER) & Pirrong (2012,
    // Commodity Price Dynamics: A Structural Approach, Ch. 2-3) Working Storage Supply Curve,
    // Stockout Convenience Yield Hyperbola & Contango-to-Backwardation Regime Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct WorkingBrennanStorageConvenienceYieldSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub inventory_to_consumption_ratio_i: T,
        pub minimum_pipeline_inventory_i_min: T,
        pub convenience_yield_convexity_alpha: T,
        pub stockout_steepness_exponent_gamma: T,
        pub physical_storage_cost_rate_u: T,
        pub risk_free_rate_r: T,
        pub futures_maturity_years_t: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct WorkingBrennanStorageConvenienceYieldReport<T: Real = f64> {
        pub implied_convenience_yield_y: T,
        pub net_cost_of_carry_rate_b: T,
        pub working_marginal_price_of_storage_slope: T,
        pub equilibrium_futures_price: T,
        pub annualized_futures_basis_pct: T,
        pub critical_backwardation_inventory_threshold: T,
    }

    pub fn evaluate_working_brennan_storage_convenience_yield_curve<T: Real>(
        spec: &WorkingBrennanStorageConvenienceYieldSpec<T>,
    ) -> WorkingBrennanStorageConvenienceYieldReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1.0));
        let inv = spec
            .inventory_to_consumption_ratio_i
            .max_val(T::from_f64(0.01));
        let inv_min = spec
            .minimum_pipeline_inventory_i_min
            .max_val(T::from_f64(0.005));
        let alpha = spec
            .convenience_yield_convexity_alpha
            .max_val(T::from_f64(1e-4));
        let gamma = spec
            .stockout_steepness_exponent_gamma
            .max_val(T::from_f64(0.10));
        let u = spec.physical_storage_cost_rate_u.max_val(T::zero());
        let r = spec.risk_free_rate_r;
        let t = spec.futures_maturity_years_t.max_val(T::from_f64(0.05));

        let diff = inv - inv_min;
        let rad = (diff * diff + T::from_f64(1e-4)).sqrt();
        let excess_buf = T::from_f64(0.5) * (rad + diff) + T::from_f64(0.02);

        let y_conv = alpha * excess_buf.pow(-gamma);
        let dbuf_dinv = T::from_f64(0.5) * (T::one() + diff / rad);
        let working_slope = alpha * gamma * excess_buf.pow(-gamma - T::one()) * dbuf_dinv;

        let b_carry = r + u - y_conv;
        let f_eq = s0 * (b_carry * t).exp();
        let basis_pct = ((f_eq - s0) / (s0 * t)) * T::from_f64(100.0);

        let carry_gross = (r + u).max_val(T::from_f64(1e-4));
        let inv_crit = inv_min + (alpha / carry_gross).pow(T::one() / gamma) - T::from_f64(0.02);

        WorkingBrennanStorageConvenienceYieldReport {
            implied_convenience_yield_y: y_conv,
            net_cost_of_carry_rate_b: b_carry,
            working_marginal_price_of_storage_slope: working_slope,
            equilibrium_futures_price: f_eq,
            annualized_futures_basis_pct: basis_pct,
            critical_backwardation_inventory_threshold: inv_crit,
        }
    }

    #[cfg(test)]
    mod evaluate_working_brennan_storage_convenience_yield_curve_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn working_brennan_storage_convenience_yield_invariants_and_ad_parity() {
            let spec = WorkingBrennanStorageConvenienceYieldSpec {
                spot_price_s0: 75.0,
                inventory_to_consumption_ratio_i: 0.11,
                minimum_pipeline_inventory_i_min: 0.08,
                convenience_yield_convexity_alpha: 0.018,
                stockout_steepness_exponent_gamma: 0.85,
                physical_storage_cost_rate_u: 0.045,
                risk_free_rate_r: 0.040,
                futures_maturity_years_t: 1.0,
            };
            let rep = evaluate_working_brennan_storage_convenience_yield_curve(&spec);
            assert!(rep.implied_convenience_yield_y > 0.10);
            assert!(rep.net_cost_of_carry_rate_b < 0.0);
            assert!(rep.equilibrium_futures_price < spec.spot_price_s0);
            assert!(rep.working_marginal_price_of_storage_slope > 0.0);

            let h = 1e-5;
            let r_ad = evaluate_working_brennan_storage_convenience_yield_curve(
                &WorkingBrennanStorageConvenienceYieldSpec {
                    spot_price_s0: Dual::<8>::variable(spec.spot_price_s0, 0),
                    inventory_to_consumption_ratio_i: Dual::<8>::variable(
                        spec.inventory_to_consumption_ratio_i,
                        1,
                    ),
                    minimum_pipeline_inventory_i_min: Dual::<8>::variable(
                        spec.minimum_pipeline_inventory_i_min,
                        2,
                    ),
                    convenience_yield_convexity_alpha: Dual::<8>::variable(
                        spec.convenience_yield_convexity_alpha,
                        3,
                    ),
                    stockout_steepness_exponent_gamma: Dual::<8>::variable(
                        spec.stockout_steepness_exponent_gamma,
                        4,
                    ),
                    physical_storage_cost_rate_u: Dual::<8>::variable(
                        spec.physical_storage_cost_rate_u,
                        5,
                    ),
                    risk_free_rate_r: Dual::<8>::variable(spec.risk_free_rate_r, 6),
                    futures_maturity_years_t: Dual::<8>::variable(spec.futures_maturity_years_t, 7),
                },
            );
            for (idx, &ad) in r_ad.equilibrium_futures_price.d.iter().enumerate() {
                let mut sp = spec;
                let mut sm = spec;
                match idx {
                    0 => {
                        sp.spot_price_s0 += h;
                        sm.spot_price_s0 -= h;
                    }
                    1 => {
                        sp.inventory_to_consumption_ratio_i += h;
                        sm.inventory_to_consumption_ratio_i -= h;
                    }
                    2 => {
                        sp.minimum_pipeline_inventory_i_min += h;
                        sm.minimum_pipeline_inventory_i_min -= h;
                    }
                    3 => {
                        sp.convenience_yield_convexity_alpha += h;
                        sm.convenience_yield_convexity_alpha -= h;
                    }
                    4 => {
                        sp.stockout_steepness_exponent_gamma += h;
                        sm.stockout_steepness_exponent_gamma -= h;
                    }
                    5 => {
                        sp.physical_storage_cost_rate_u += h;
                        sm.physical_storage_cost_rate_u -= h;
                    }
                    6 => {
                        sp.risk_free_rate_r += h;
                        sm.risk_free_rate_r -= h;
                    }
                    _ => {
                        sp.futures_maturity_years_t += h;
                        sm.futures_maturity_years_t -= h;
                    }
                }
                let fd = (evaluate_working_brennan_storage_convenience_yield_curve(&sp)
                    .equilibrium_futures_price
                    - evaluate_working_brennan_storage_convenience_yield_curve(&sm)
                        .equilibrium_futures_price)
                    / (2.0 * h);
                assert!(ad.abs() > 1e-6);
                assert!((ad - fd).abs() < 1e-3);
            }
        }
    }
}
pub use evaluate_working_brennan_storage_convenience_yield_curve_module::*;
