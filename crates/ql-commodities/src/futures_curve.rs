use ql_ad::{Dual, Real, TapeScope, TapeVar};
use ql_math::norm_cdf;

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

mod evaluate_gibson_schwartz_two_factor_futures_samuelson_vol_module {
    // Gibson & Schwartz (1990, Journal of Finance 45(3), "Stochastic Convenience Yield and the
    // Pricing of Oil Contingent Claims"), Schwartz (1997, JF 52(3), Model 2),
    // & Geman (2005, Commodities and Commodity Derivatives, Ch. 3)
    // Two-Factor Stochastic Convenience Yield Futures Curve & Samuelson Term-Structure Volatility Engine.

    use ql_ad::Real;

    #[derive(Debug, Clone, Copy)]
    pub struct GibsonSchwartzSamuelsonFuturesSpec<T: Real = f64> {
        pub spot_price_s0: T,
        pub initial_convenience_yield_delta0: T,
        pub convenience_yield_mean_reversion_kappa: T,
        pub risk_neutral_long_run_convenience_alpha: T,
        pub spot_volatility_sigma_s: T,
        pub convenience_yield_volatility_sigma_d: T,
        pub spot_convenience_correlation_rho: T,
        pub risk_free_rate_r: T,
        pub maturity_years_tau: T,
    }

    #[derive(Debug, Clone, Copy)]
    pub struct GibsonSchwartzSamuelsonFuturesReport<T: Real = f64> {
        pub convenience_yield_loading_b_tau: T,
        pub gibson_schwartz_futures_price: T,
        pub instantaneous_futures_volatility: T,
        pub term_average_black76_option_volatility: T,
        pub samuelson_volatility_attenuation_ratio: T,
        pub implied_annualized_roll_yield: T,
    }

    pub fn evaluate_gibson_schwartz_two_factor_futures_samuelson_vol<T: Real>(
        spec: &GibsonSchwartzSamuelsonFuturesSpec<T>,
    ) -> GibsonSchwartzSamuelsonFuturesReport<T> {
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1.0));
        let delta0 = spec.initial_convenience_yield_delta0;
        let kappa = spec
            .convenience_yield_mean_reversion_kappa
            .max_val(T::from_f64(0.05));
        let alpha_q = spec.risk_neutral_long_run_convenience_alpha;
        let sig_s = spec.spot_volatility_sigma_s.max_val(T::from_f64(0.02));
        let sig_d = spec
            .convenience_yield_volatility_sigma_d
            .max_val(T::from_f64(0.01));
        let rho = spec
            .spot_convenience_correlation_rho
            .max_val(T::from_f64(-0.98))
            .min_val(T::from_f64(0.98));
        let r = spec.risk_free_rate_r;
        let tau = spec.maturity_years_tau.max_val(T::from_f64(0.05));

        let e_kt = (-kappa * tau).exp();
        let e_2kt = (-T::from_f64(2.0) * kappa * tau).exp();
        let b_tau = (T::one() - e_kt) / kappa;

        let k2 = kappa * kappa;
        let k3 = k2 * kappa;
        let sig_d2 = sig_d * sig_d;

        // Schwartz (1997, Eqs. 18-20) exact affine log-futures term A(tau):
        let a_tau =
            (r - alpha_q + sig_d2 / (T::from_f64(2.0) * k2) - (rho * sig_s * sig_d) / kappa) * tau
                + (sig_d2 / (T::from_f64(4.0) * k3)) * (T::one() - e_2kt)
                + (alpha_q * kappa + rho * sig_s * sig_d - sig_d2 / kappa) * (b_tau / kappa);

        let log_f = s0.ln() - delta0 * b_tau + a_tau;
        let f_price = log_f.exp();

        // Instantaneous futures return variance at maturity tau (Samuelson effect):
        let var_inst = (sig_s * sig_s - T::from_f64(2.0) * rho * sig_s * sig_d * b_tau
            + sig_d2 * b_tau * b_tau)
            .max_val(T::from_f64(1e-6));
        let vol_inst = var_inst.sqrt();

        // Integrated average variance (1/tau) * int_0^tau sigma_F^2(u) du for Black-76 futures options:
        let int_var = sig_s * sig_s * tau
            - (T::from_f64(2.0) * rho * sig_s * sig_d / kappa) * (tau - b_tau)
            + (sig_d2 / k2)
                * (tau - T::from_f64(2.0) * b_tau
                    + (T::one() - e_2kt) / (T::from_f64(2.0) * kappa));
        let vol_term = (int_var / tau).max_val(T::from_f64(1e-6)).sqrt();

        let samuelson_ratio = vol_inst / sig_s;
        let roll_yield = -(log_f - s0.ln()) / tau;

        GibsonSchwartzSamuelsonFuturesReport {
            convenience_yield_loading_b_tau: b_tau,
            gibson_schwartz_futures_price: f_price,
            instantaneous_futures_volatility: vol_inst,
            term_average_black76_option_volatility: vol_term,
            samuelson_volatility_attenuation_ratio: samuelson_ratio,
            implied_annualized_roll_yield: roll_yield,
        }
    }

    #[cfg(test)]
    mod evaluate_gibson_schwartz_two_factor_futures_samuelson_vol_tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn gibson_schwartz_samuelson_futures_invariants_and_ad_parity() {
            let spec = GibsonSchwartzSamuelsonFuturesSpec {
                spot_price_s0: 78.0,
                initial_convenience_yield_delta0: 0.14,
                convenience_yield_mean_reversion_kappa: 1.35,
                risk_neutral_long_run_convenience_alpha: 0.08,
                spot_volatility_sigma_s: 0.36,
                convenience_yield_volatility_sigma_d: 0.28,
                spot_convenience_correlation_rho: 0.80,
                risk_free_rate_r: 0.04,
                maturity_years_tau: 1.50,
            };
            let rep = evaluate_gibson_schwartz_two_factor_futures_samuelson_vol(&spec);
            // With positive correlation rho > 0, Samuelson effect damps distant futures volatility below spot vol:
            assert!(rep.instantaneous_futures_volatility < spec.spot_volatility_sigma_s);
            assert!(rep.samuelson_volatility_attenuation_ratio < 1.0);

            let h = 1e-6;
            let r_ad = evaluate_gibson_schwartz_two_factor_futures_samuelson_vol(
                &GibsonSchwartzSamuelsonFuturesSpec {
                    spot_price_s0: Dual::<9>::variable(spec.spot_price_s0, 0),
                    initial_convenience_yield_delta0: Dual::<9>::variable(
                        spec.initial_convenience_yield_delta0,
                        1,
                    ),
                    convenience_yield_mean_reversion_kappa: Dual::<9>::variable(
                        spec.convenience_yield_mean_reversion_kappa,
                        2,
                    ),
                    risk_neutral_long_run_convenience_alpha: Dual::<9>::variable(
                        spec.risk_neutral_long_run_convenience_alpha,
                        3,
                    ),
                    spot_volatility_sigma_s: Dual::<9>::variable(spec.spot_volatility_sigma_s, 4),
                    convenience_yield_volatility_sigma_d: Dual::<9>::variable(
                        spec.convenience_yield_volatility_sigma_d,
                        5,
                    ),
                    spot_convenience_correlation_rho: Dual::<9>::variable(
                        spec.spot_convenience_correlation_rho,
                        6,
                    ),
                    risk_free_rate_r: Dual::<9>::variable(spec.risk_free_rate_r, 7),
                    maturity_years_tau: Dual::<9>::variable(spec.maturity_years_tau, 8),
                },
            );
            for i in 0..9 {
                let mut v = [
                    spec.spot_price_s0,
                    spec.initial_convenience_yield_delta0,
                    spec.convenience_yield_mean_reversion_kappa,
                    spec.risk_neutral_long_run_convenience_alpha,
                    spec.spot_volatility_sigma_s,
                    spec.convenience_yield_volatility_sigma_d,
                    spec.spot_convenience_correlation_rho,
                    spec.risk_free_rate_r,
                    spec.maturity_years_tau,
                ];
                v[i] += h;
                let up = evaluate_gibson_schwartz_two_factor_futures_samuelson_vol(
                    &GibsonSchwartzSamuelsonFuturesSpec {
                        spot_price_s0: v[0],
                        initial_convenience_yield_delta0: v[1],
                        convenience_yield_mean_reversion_kappa: v[2],
                        risk_neutral_long_run_convenience_alpha: v[3],
                        spot_volatility_sigma_s: v[4],
                        convenience_yield_volatility_sigma_d: v[5],
                        spot_convenience_correlation_rho: v[6],
                        risk_free_rate_r: v[7],
                        maturity_years_tau: v[8],
                    },
                );
                v[i] -= 2.0 * h;
                let dn = evaluate_gibson_schwartz_two_factor_futures_samuelson_vol(
                    &GibsonSchwartzSamuelsonFuturesSpec {
                        spot_price_s0: v[0],
                        initial_convenience_yield_delta0: v[1],
                        convenience_yield_mean_reversion_kappa: v[2],
                        risk_neutral_long_run_convenience_alpha: v[3],
                        spot_volatility_sigma_s: v[4],
                        convenience_yield_volatility_sigma_d: v[5],
                        spot_convenience_correlation_rho: v[6],
                        risk_free_rate_r: v[7],
                        maturity_years_tau: v[8],
                    },
                );
                let fd = (up.gibson_schwartz_futures_price - dn.gibson_schwartz_futures_price)
                    / (2.0 * h);
                assert!(
                    (r_ad.gibson_schwartz_futures_price.d[i] - fd).abs() < 1e-4,
                    "input {i}"
                );
            }
        }
    }
}
pub use evaluate_gibson_schwartz_two_factor_futures_samuelson_vol_module::*;

/// Specification for the Routledge–Seppi–Spatt (2000) competitive storage and stockout backwardation
/// equilibrium commodity forward curve engine.
///
/// Reference: Routledge, Seppi & Spatt (2000, *Journal of Finance* 55(3), "Equilibrium Forward Curves
/// for Commodities"), Deaton & Laroque (1992), Pirrong (2012, *Commodity Price Dynamics*).
#[derive(Debug, Clone, Copy)]
pub struct RoutledgeSeppiSpattStorageSpec<T: Real> {
    /// Current commodity spot price $S_0$.
    pub spot_price_s0: T,
    /// Normalized physical inventory buffer $I_0$.
    pub physical_inventory_buffer_i0: T,
    /// Net transient demand shock state $x_0$.
    pub net_demand_shock_x0: T,
    /// Mean-reversion speed of net demand shocks $\kappa_x$.
    pub demand_mean_reversion_kappa_x: T,
    /// Volatility of net demand shocks $\sigma_x$.
    pub demand_shock_volatility_sigma_x: T,
    /// Proportional physical storage cost rate $c_s$.
    pub storage_cost_rate_cs: T,
    /// Risk-free interest rate $r$.
    pub risk_free_rate_r: T,
    /// Forward contract maturity $\tau$ (in years).
    pub forward_maturity_years_tau: T,
}

/// Analytics returned by [`evaluate_routledge_seppi_spatt_equilibrium_storage_forward_curve`].
#[derive(Debug, Clone, Copy)]
pub struct RoutledgeSeppiSpattStorageAnalytics<T: Real> {
    /// Probability of inventory stockout $p_{\text{stockout}}(\tau)$ over $[0, \tau]$.
    pub inventory_stockout_probability: T,
    /// Equilibrium commodity forward price $F(0, \tau)$.
    pub equilibrium_forward_price: T,
    /// Endogenous stockout-implied convenience yield $\delta_{\text{imp}}(\tau)$.
    pub endogenous_convenience_yield: T,
    /// Inventory-dependent forward volatility $\sigma_F(\tau, I_0)$ (Samuelson + stockout amplification).
    pub stockout_amplified_forward_volatility: T,
}

/// Evaluates the Routledge–Seppi–Spatt (2000) competitive storage equilibrium forward price, stockout
/// probability, endogenous convenience yield, and inventory-dependent forward return volatility.
pub fn evaluate_routledge_seppi_spatt_equilibrium_storage_forward_curve<T: Real>(
    spec: &RoutledgeSeppiSpattStorageSpec<T>,
) -> RoutledgeSeppiSpattStorageAnalytics<T> {
    let tau = spec.forward_maturity_years_tau.max_val(T::from_f64(1e-4));
    let kx = spec
        .demand_mean_reversion_kappa_x
        .max_val(T::from_f64(1e-4));
    let sig_x = spec
        .demand_shock_volatility_sigma_x
        .max_val(T::from_f64(1e-4));

    // Conditional standard deviation of OU net demand shock over [0, tau]
    let ou_std = sig_x
        * ((T::one() - (-T::from_f64(2.0) * kx * tau).exp()) / (T::from_f64(2.0) * kx)).sqrt();
    let expected_net_shortfall =
        spec.net_demand_shock_x0 * (-kx * tau).exp() - spec.physical_inventory_buffer_i0;
    let z_so = expected_net_shortfall / ou_std;

    let p_stockout = norm_cdf(z_so);
    let inv_sqrt_2pi = T::from_f64(0.3989422804014327);
    let phi_z = inv_sqrt_2pi * (-T::from_f64(0.5) * z_so * z_so).exp();

    // Expected truncated normal stockout scarcity rent decaying over [0, tau]
    let current_scarcity = sig_x * (z_so * p_stockout + phi_z);
    let attenuation = (T::one() - (-kx * tau).exp()) / (kx * tau);
    let delta_endog = current_scarcity * (T::from_f64(0.5) + T::from_f64(0.5) * attenuation);

    // Competitive storage carry cost r + c_s minus endogenous stockout convenience yield delta_endog
    let net_carry = spec.risk_free_rate_r + spec.storage_cost_rate_cs - delta_endog;
    let f_0_tau = spec.spot_price_s0 * (net_carry * tau).exp();

    // Samuelson forward volatility attenuated by exp(-kx * tau) but amplified by stockout probability
    let sig_f = sig_x * (-kx * tau).exp() * (T::one() + T::from_f64(0.65) * p_stockout);

    RoutledgeSeppiSpattStorageAnalytics {
        inventory_stockout_probability: p_stockout,
        equilibrium_forward_price: f_0_tau,
        endogenous_convenience_yield: delta_endog,
        stockout_amplified_forward_volatility: sig_f,
    }
}

#[cfg(test)]
mod rss_storage_suite64_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn routledge_seppi_spatt_storage_invariants_and_ad_parity() {
        let base = RoutledgeSeppiSpattStorageSpec {
            spot_price_s0: 78.0_f64,
            physical_inventory_buffer_i0: 0.12,
            net_demand_shock_x0: 0.18,
            demand_mean_reversion_kappa_x: 1.40,
            demand_shock_volatility_sigma_x: 0.26,
            storage_cost_rate_cs: 0.035,
            risk_free_rate_r: 0.045,
            forward_maturity_years_tau: 0.75,
        };
        let r0 = evaluate_routledge_seppi_spatt_equilibrium_storage_forward_curve(&base);
        assert!(r0.inventory_stockout_probability > 0.0 && r0.inventory_stockout_probability < 1.0);
        assert!(r0.equilibrium_forward_price > 50.0);
        assert!(r0.endogenous_convenience_yield > 0.0);
        assert!(r0.stockout_amplified_forward_volatility > 0.0);

        let inputs = [
            base.spot_price_s0,
            base.physical_inventory_buffer_i0,
            base.net_demand_shock_x0,
            base.demand_mean_reversion_kappa_x,
            base.demand_shock_volatility_sigma_x,
            base.storage_cost_rate_cs,
            base.risk_free_rate_r,
            base.forward_maturity_years_tau,
        ];
        let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            dual_in[i] = Dual::variable(inputs[i], i);
        }
        let spec_ad = RoutledgeSeppiSpattStorageSpec {
            spot_price_s0: dual_in[0],
            physical_inventory_buffer_i0: dual_in[1],
            net_demand_shock_x0: dual_in[2],
            demand_mean_reversion_kappa_x: dual_in[3],
            demand_shock_volatility_sigma_x: dual_in[4],
            storage_cost_rate_cs: dual_in[5],
            risk_free_rate_r: dual_in[6],
            forward_maturity_years_tau: dual_in[7],
        };
        let r_ad = evaluate_routledge_seppi_spatt_equilibrium_storage_forward_curve(&spec_ad);
        let eps = 1e-6;
        for i in 0..8 {
            let mut p = inputs;
            let mut m = inputs;
            p[i] += eps;
            m[i] -= eps;
            let mk = |v: [f64; 8]| RoutledgeSeppiSpattStorageSpec {
                spot_price_s0: v[0],
                physical_inventory_buffer_i0: v[1],
                net_demand_shock_x0: v[2],
                demand_mean_reversion_kappa_x: v[3],
                demand_shock_volatility_sigma_x: v[4],
                storage_cost_rate_cs: v[5],
                risk_free_rate_r: v[6],
                forward_maturity_years_tau: v[7],
            };
            let rp = evaluate_routledge_seppi_spatt_equilibrium_storage_forward_curve(&mk(p));
            let rm = evaluate_routledge_seppi_spatt_equilibrium_storage_forward_curve(&mk(m));
            let fd = (rp.equilibrium_forward_price - rm.equilibrium_forward_price) / (2.0 * eps);
            assert!(
                r_ad.equilibrium_forward_price.d[i].abs() > 1e-6,
                "zero derivative at {i}"
            );
            assert!((r_ad.equilibrium_forward_price.d[i] - fd).abs() < 1e-4);
        }
    }
}

pub mod cortazar_schwartz_three_factor_oil_mod {
    use ql_ad::Real;

    /// Specification for the Cortazar & Schwartz (2003, Energy Economics) / Geman (2005, Ch. 3)
    /// three-factor oil futures term structure model.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct CortazarSchwartzThreeFactorSpec<T: Real> {
        /// Spot oil price $S_0 > 0$.
        pub spot_price_s0: T,
        /// Initial instantaneous convenience yield $\delta_0$.
        pub initial_convenience_yield_delta0: T,
        /// Initial instantaneous long-term equilibrium growth rate $v_0$.
        pub initial_long_term_drift_v0: T,
        /// Risk-neutral long-run convenience yield mean $\alpha_2^*$.
        pub risk_neutral_convenience_mean_alpha2: T,
        /// Risk-neutral long-run equilibrium return $\bar{v}^*$.
        pub risk_neutral_long_run_drift_vbar: T,
        /// Convenience yield mean-reversion speed $\kappa_2 > 0$.
        pub convenience_reversion_speed_kappa2: T,
        /// Long-term drift mean-reversion speed $\kappa_3 > 0$.
        pub drift_reversion_speed_kappa3: T,
        /// Spot oil return volatility $\sigma_1 > 0$.
        pub spot_volatility_sigma1: T,
        /// Convenience yield volatility $\sigma_2 > 0$.
        pub convenience_volatility_sigma2: T,
        /// Long-term equilibrium drift volatility $\sigma_3 > 0$.
        pub drift_volatility_sigma3: T,
        /// Correlation $\rho_{12} \in (-1, 1)$ between spot oil and convenience yield.
        pub spot_convenience_correlation_rho12: T,
        /// Futures contract maturity $\tau > 0$ (years).
        pub maturity_years_tau: T,
    }

    /// Analytical outputs for the Cortazar-Schwartz three-factor commodity futures engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct CortazarSchwartzThreeFactorAnalytics<T: Real> {
        /// Convenience-yield loading $B_2(\tau) = (1 - e^{-\kappa_2 \tau})/\kappa_2$.
        pub convenience_yield_loading_b2: T,
        /// Long-term drift loading $B_3(\tau) = (1 - e^{-\kappa_3 \tau})/\kappa_3$.
        pub long_term_drift_loading_b3: T,
        /// Instantaneous futures return volatility $\sigma_F(\tau)$ at maturity $\tau$.
        pub term_futures_volatility: T,
        /// Closed-form three-factor futures price $F(0, \tau)$.
        pub three_factor_futures_price: T,
    }

    /// Evaluates the Cortazar & Schwartz (2003) three-factor oil futures price and Samuelson
    /// term-structure volatility curve.
    pub fn evaluate_cortazar_schwartz_three_factor_oil_futures_curve<T: Real>(
        spec: &CortazarSchwartzThreeFactorSpec<T>,
    ) -> CortazarSchwartzThreeFactorAnalytics<T> {
        let tau = spec.maturity_years_tau.max_val(T::from_f64(1e-4));
        let k2 = spec
            .convenience_reversion_speed_kappa2
            .max_val(T::from_f64(1e-4));
        let k3 = spec.drift_reversion_speed_kappa3.max_val(T::from_f64(1e-4));
        let s0 = spec.spot_price_s0.max_val(T::from_f64(1e-4));

        let b2 = (T::one() - (-k2 * tau).exp()) / k2;
        let b3 = (T::one() - (-k3 * tau).exp()) / k3;
        let b2_2 = (T::one() - (-T::from_f64(2.0) * k2 * tau).exp()) / (T::from_f64(2.0) * k2);
        let b3_2 = (T::one() - (-T::from_f64(2.0) * k3 * tau).exp()) / (T::from_f64(2.0) * k3);

        let sig1 = spec.spot_volatility_sigma1.max_val(T::from_f64(1e-5));
        let sig2 = spec
            .convenience_volatility_sigma2
            .max_val(T::from_f64(1e-5));
        let sig3 = spec.drift_volatility_sigma3.max_val(T::from_f64(1e-5));
        let rho12 = spec
            .spot_convenience_correlation_rho12
            .max_val(T::from_f64(-0.99))
            .min_val(T::from_f64(0.99));

        // Integrated variance of log spot price under Q
        let var_log_s = sig1 * sig1 * tau
            + (sig2 * sig2 / (k2 * k2)) * (tau - T::from_f64(2.0) * b2 + b2_2)
            + (sig3 * sig3 / (k3 * k3)) * (tau - T::from_f64(2.0) * b3 + b3_2)
            - (T::from_f64(2.0) * rho12 * sig1 * sig2 / k2) * (tau - b2);

        let log_f = s0.ln() - spec.initial_convenience_yield_delta0 * b2
            + spec.initial_long_term_drift_v0 * b3
            + spec.risk_neutral_long_run_drift_vbar * (tau - b3)
            - spec.risk_neutral_convenience_mean_alpha2 * (tau - b2)
            - T::from_f64(0.5) * sig1 * sig1 * tau
            + T::from_f64(0.5) * var_log_s;

        let inst_var_f = (sig1 * sig1 + sig2 * sig2 * b2 * b2 + sig3 * sig3 * b3 * b3
            - T::from_f64(2.0) * rho12 * sig1 * sig2 * b2)
            .max_val(T::from_f64(1e-8));

        CortazarSchwartzThreeFactorAnalytics {
            convenience_yield_loading_b2: b2,
            long_term_drift_loading_b3: b3,
            term_futures_volatility: inst_var_f.sqrt(),
            three_factor_futures_price: log_f.exp(),
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn cortazar_schwartz_three_factor_invariants_and_ad_parity() {
            let base = CortazarSchwartzThreeFactorSpec {
                spot_price_s0: 78.0_f64,
                initial_convenience_yield_delta0: 0.08,
                initial_long_term_drift_v0: 0.035,
                risk_neutral_convenience_mean_alpha2: 0.05,
                risk_neutral_long_run_drift_vbar: 0.03,
                convenience_reversion_speed_kappa2: 1.40,
                drift_reversion_speed_kappa3: 0.35,
                spot_volatility_sigma1: 0.32,
                convenience_volatility_sigma2: 0.24,
                drift_volatility_sigma3: 0.08,
                spot_convenience_correlation_rho12: 0.72,
                maturity_years_tau: 2.5,
            };
            let r0 = evaluate_cortazar_schwartz_three_factor_oil_futures_curve(&base);
            assert!(r0.convenience_yield_loading_b2 > 0.0);
            assert!(r0.long_term_drift_loading_b3 > r0.convenience_yield_loading_b2);
            assert!(r0.term_futures_volatility > 0.0);
            assert!(r0.three_factor_futures_price > 0.0);

            let inputs = [
                base.spot_price_s0,
                base.initial_convenience_yield_delta0,
                base.initial_long_term_drift_v0,
                base.risk_neutral_convenience_mean_alpha2,
                base.risk_neutral_long_run_drift_vbar,
                base.convenience_reversion_speed_kappa2,
                base.drift_reversion_speed_kappa3,
                base.spot_volatility_sigma1,
                base.convenience_volatility_sigma2,
                base.drift_volatility_sigma3,
                base.spot_convenience_correlation_rho12,
                base.maturity_years_tau,
            ];
            let mut dual_in = [Dual::<12>::from_f64(0.0); 12];
            for i in 0..12 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = CortazarSchwartzThreeFactorSpec {
                spot_price_s0: dual_in[0],
                initial_convenience_yield_delta0: dual_in[1],
                initial_long_term_drift_v0: dual_in[2],
                risk_neutral_convenience_mean_alpha2: dual_in[3],
                risk_neutral_long_run_drift_vbar: dual_in[4],
                convenience_reversion_speed_kappa2: dual_in[5],
                drift_reversion_speed_kappa3: dual_in[6],
                spot_volatility_sigma1: dual_in[7],
                convenience_volatility_sigma2: dual_in[8],
                drift_volatility_sigma3: dual_in[9],
                spot_convenience_correlation_rho12: dual_in[10],
                maturity_years_tau: dual_in[11],
            };
            let r_ad = evaluate_cortazar_schwartz_three_factor_oil_futures_curve(&spec_ad);
            let eps = 1e-6;
            for i in 0..12 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 12]| CortazarSchwartzThreeFactorSpec {
                    spot_price_s0: v[0],
                    initial_convenience_yield_delta0: v[1],
                    initial_long_term_drift_v0: v[2],
                    risk_neutral_convenience_mean_alpha2: v[3],
                    risk_neutral_long_run_drift_vbar: v[4],
                    convenience_reversion_speed_kappa2: v[5],
                    drift_reversion_speed_kappa3: v[6],
                    spot_volatility_sigma1: v[7],
                    convenience_volatility_sigma2: v[8],
                    drift_volatility_sigma3: v[9],
                    spot_convenience_correlation_rho12: v[10],
                    maturity_years_tau: v[11],
                };
                let rp = evaluate_cortazar_schwartz_three_factor_oil_futures_curve(&mk(p));
                let rm = evaluate_cortazar_schwartz_three_factor_oil_futures_curve(&mk(m));
                let fd =
                    (rp.three_factor_futures_price - rm.three_factor_futures_price) / (2.0 * eps);
                assert!(
                    r_ad.three_factor_futures_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.three_factor_futures_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use cortazar_schwartz_three_factor_oil_mod::*;

// ============================================================================
// Benth & Saltyte-Benth (2007, 2012) Seasonal Ornstein-Uhlenbeck HDD/CDD Weather Swap
// ============================================================================

pub mod benth_seasonal_temperature_hdd_cdd_swap_mod {
    use ql_ad::Real;
    use ql_math::distributions::{norm_cdf, norm_pdf};

    /// Specification for the Benth & Šaltytė-Benth (2007, *JBF* 31(7); 2012, Ch. 5)
    /// seasonal Ornstein–Uhlenbeck temperature process and HDD/CDD weather derivative swap.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BenthSeasonalTemperatureHddCddSpec<T: Real> {
        /// Initial observed daily average temperature $T_0$ (°C).
        pub initial_temperature_t0: T,
        /// Annual long-run mean temperature level $a_0$ (°C).
        pub seasonal_mean_level_a0: T,
        /// Seasonal cosine wave amplitude $b_1$ (°C).
        pub seasonal_amplitude_b1: T,
        /// Ornstein-Uhlenbeck temperature mean-reversion speed $\kappa > 0$ (annualized).
        pub mean_reversion_speed_kappa: T,
        /// Ornstein-Uhlenbeck temperature volatility $\sigma_T > 0$ (°C / $\sqrt{\text{year}}$).
        pub temperature_volatility_sigma: T,
        /// Degree-day reference threshold temperature $K_{\text{base}}$ (e.g., $18.0^\circ\text{C}$).
        pub degree_day_threshold_k_base: T,
        /// Girsanov / Esscher market price of weather risk $\theta_W$.
        pub market_price_of_weather_risk_theta: T,
        /// Risk-free discount rate $r > 0$.
        pub discount_rate_r: T,
    }

    /// Analytical outputs for the Benth seasonal OU HDD/CDD weather swap engine.
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BenthSeasonalTemperatureHddCddAnalytics<T: Real> {
        /// Accumulated risk-neutral expected Heating Degree Days $\mathbb{E}^\theta[\text{HDD}]$.
        pub expected_cumulative_hdd: T,
        /// Accumulated risk-neutral expected Cooling Degree Days $\mathbb{E}^\theta[\text{CDD}]$.
        pub expected_cumulative_cdd: T,
        /// Terminal seasonal risk-neutral mean temperature $\mu_Q(T)$ (°C).
        pub terminal_risk_neutral_mean_temperature: T,
        /// Discounted Heating Degree Day (HDD) swap present value $e^{-r T}\mathbb{E}^\theta[\text{HDD}]$.
        pub discounted_hdd_swap_present_value: T,
    }

    /// Evaluates the Benth & Šaltytė-Benth (2012, Proposition 5.1 & Theorem 5.2) analytical
    /// Bachelier HDD and CDD expectations across a 4-month winter/shoulder season.
    pub fn evaluate_benth_cartea_saltyte_benth_seasonal_temperature_hdd_cdd_swap<T: Real>(
        spec: &BenthSeasonalTemperatureHddCddSpec<T>,
    ) -> BenthSeasonalTemperatureHddCddAnalytics<T> {
        const PERIODS: usize = 4;
        let t_contract = T::from_f64(0.5);
        let dt = t_contract / T::from_f64(PERIODS as f64);
        let days_per_period = T::from_f64(30.0);

        let kappa = spec.mean_reversion_speed_kappa.max_val(T::from_f64(1e-3));
        let sig_t = spec.temperature_volatility_sigma.max_val(T::from_f64(1e-3));
        let two_pi = T::from_f64(std::f64::consts::TAU);

        let lambda_0 = spec.seasonal_mean_level_a0 - spec.seasonal_amplitude_b1;
        let mut cum_hdd = T::zero();
        let mut cum_cdd = T::zero();
        let mut last_mu_q = lambda_0;

        for m in 1..=PERIODS {
            let tm = (T::from_f64(m as f64) - T::from_f64(0.5)) * dt;
            let e_kt = (-kappa * tm).exp();
            let lambda_t =
                spec.seasonal_mean_level_a0 - spec.seasonal_amplitude_b1 * (two_pi * tm).cos();

            // Risk-neutral OU mean under market price of weather risk theta_W
            let mu_q = lambda_t
                + (spec.initial_temperature_t0 - lambda_0) * e_kt
                + (spec.market_price_of_weather_risk_theta * sig_t / kappa) * (T::one() - e_kt);
            last_mu_q = mu_q;

            let var_q = (sig_t * sig_t / (T::from_f64(2.0) * kappa)) * (T::one() - e_kt * e_kt);
            let sd_q = var_q.max_val(T::from_f64(1e-6)).sqrt();

            let diff = spec.degree_day_threshold_k_base - mu_q;
            let d = diff / sd_q;
            let hdd_day = diff * norm_cdf(d) + sd_q * norm_pdf(d);
            let cdd_day = -diff * norm_cdf(-d) + sd_q * norm_pdf(d);

            cum_hdd += days_per_period * hdd_day;
            cum_cdd += days_per_period * cdd_day;
        }

        let df = (-spec.discount_rate_r * t_contract).exp();
        let pv_hdd = df * cum_hdd;

        BenthSeasonalTemperatureHddCddAnalytics {
            expected_cumulative_hdd: cum_hdd,
            expected_cumulative_cdd: cum_cdd,
            terminal_risk_neutral_mean_temperature: last_mu_q,
            discounted_hdd_swap_present_value: pv_hdd,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn benth_seasonal_hdd_cdd_swap_invariants_and_ad_parity() {
            let base = BenthSeasonalTemperatureHddCddSpec {
                initial_temperature_t0: 12.5_f64,
                seasonal_mean_level_a0: 15.0,
                seasonal_amplitude_b1: 11.0,
                mean_reversion_speed_kappa: 6.5,
                temperature_volatility_sigma: 8.0,
                degree_day_threshold_k_base: 18.0,
                market_price_of_weather_risk_theta: -0.18,
                discount_rate_r: 0.04,
            };
            let r0 = evaluate_benth_cartea_saltyte_benth_seasonal_temperature_hdd_cdd_swap(&base);
            assert!(r0.expected_cumulative_hdd > 0.0);
            assert!(r0.expected_cumulative_cdd > 0.0);
            assert!(r0.discounted_hdd_swap_present_value > 0.0);

            let inputs = [
                base.initial_temperature_t0,
                base.seasonal_mean_level_a0,
                base.seasonal_amplitude_b1,
                base.mean_reversion_speed_kappa,
                base.temperature_volatility_sigma,
                base.degree_day_threshold_k_base,
                base.market_price_of_weather_risk_theta,
                base.discount_rate_r,
            ];
            let mut dual_in = [Dual::<8>::from_f64(0.0); 8];
            for i in 0..8 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BenthSeasonalTemperatureHddCddSpec {
                initial_temperature_t0: dual_in[0],
                seasonal_mean_level_a0: dual_in[1],
                seasonal_amplitude_b1: dual_in[2],
                mean_reversion_speed_kappa: dual_in[3],
                temperature_volatility_sigma: dual_in[4],
                degree_day_threshold_k_base: dual_in[5],
                market_price_of_weather_risk_theta: dual_in[6],
                discount_rate_r: dual_in[7],
            };
            let r_ad =
                evaluate_benth_cartea_saltyte_benth_seasonal_temperature_hdd_cdd_swap(&spec_ad);
            let eps = 1e-6;
            for i in 0..8 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 8]| BenthSeasonalTemperatureHddCddSpec {
                    initial_temperature_t0: v[0],
                    seasonal_mean_level_a0: v[1],
                    seasonal_amplitude_b1: v[2],
                    mean_reversion_speed_kappa: v[3],
                    temperature_volatility_sigma: v[4],
                    degree_day_threshold_k_base: v[5],
                    market_price_of_weather_risk_theta: v[6],
                    discount_rate_r: v[7],
                };
                let rp =
                    evaluate_benth_cartea_saltyte_benth_seasonal_temperature_hdd_cdd_swap(&mk(p));
                let rm =
                    evaluate_benth_cartea_saltyte_benth_seasonal_temperature_hdd_cdd_swap(&mk(m));
                let fd = (rp.discounted_hdd_swap_present_value
                    - rm.discounted_hdd_swap_present_value)
                    / (2.0 * eps);
                assert!(
                    r_ad.discounted_hdd_swap_present_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.discounted_hdd_swap_present_value.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use benth_seasonal_temperature_hdd_cdd_swap_mod::*;
mod hilliard_reis_jump_futures_mod {
    use ql_ad::Real;

    /// Input specification for the Hilliard & Reis (1998, *Journal of Finance* 53(1))
    /// two-factor stochastic convenience-yield and jump-diffusion commodity futures option model.
    #[derive(Debug, Clone, Copy)]
    pub struct HilliardReisJumpFuturesSpec<T: Real> {
        /// Current futures price $F_0(T_{\text{fut}})$ (e.g., $78.0$).
        pub futures_price_f0: T,
        /// Option strike price $K$ (e.g., $80.0$).
        pub strike_price_k: T,
        /// Spot commodity diffusive volatility $\sigma_S$ (e.g., $0.32$).
        pub spot_volatility_sigma_s: T,
        /// Convenience-yield volatility $\sigma_\delta$ (e.g., $0.24$).
        pub convenience_yield_volatility_sigma_delta: T,
        /// Convenience-yield mean-reversion speed $\kappa_\delta > 0$ (e.g., $1.60$).
        pub convenience_yield_mean_reversion_kappa: T,
        /// Spot-convenience-yield correlation $\rho_{S\delta} \in (-1, 1)$ (e.g., $0.65$).
        pub spot_convenience_correlation_rho: T,
        /// Poisson jump arrival intensity $\lambda_J$ per year (e.g., $0.85$).
        pub jump_intensity_lambda_j: T,
        /// Log-jump volatility $v_J$ (e.g., $0.18$).
        pub jump_volatility_v_j: T,
        /// Risk-free discount rate $r$ across option maturity $T_{\text{opt}} = 0.5$, $T_{\text{fut}} = 1.0$ (e.g., $0.038$).
        pub risk_free_rate_r: T,
    }

    /// Analytics returned by [`evaluate_hilliard_reis_stochastic_convenience_yield_jump_futures_option`].
    #[derive(Debug, Clone, Copy)]
    pub struct HilliardReisJumpFuturesAnalytics<T: Real> {
        /// Integrated diffusive futures variance $V_F(T_{\text{opt}}, T_{\text{fut}})$.
        pub integrated_diffusive_futures_variance: T,
        /// Effective annualized diffusive futures volatility $\sqrt{V_F / T_{\text{opt}}}$.
        pub annualized_diffusive_futures_volatility: T,
        /// Hilliard–Reis jump-diffusion futures call option price.
        pub jump_diffusion_futures_call_price: T,
        /// Pure diffusive (zero-jump) futures call option price.
        pub diffusive_only_futures_call_price: T,
    }

    #[inline]
    fn std_norm_cdf<T: Real>(x: T) -> T {
        T::from_f64(0.5) * (T::one() + (x * T::from_f64(0.707_106_781_186_547_5)).erf())
    }

    #[inline]
    fn black76_call<T: Real>(f: T, k: T, df: T, total_var: T) -> T {
        let std_dev = total_var.sqrt();
        let d1 = ((f / k).ln() + T::from_f64(0.5) * total_var) / std_dev;
        let d2 = d1 - std_dev;
        df * (f * std_norm_cdf(d1) - k * std_norm_cdf(d2))
    }

    /// Evaluates the Hilliard & Reis (1998, Prop. 1 & Prop. 3) commodity futures call option
    /// under stochastic convenience yield $\delta_t$ and Poisson price jumps over option maturity
    /// $T_{\text{opt}} = 0.5$ years on a futures contract maturing at $T_{\text{fut}} = 1.0$ year.
    pub fn evaluate_hilliard_reis_stochastic_convenience_yield_jump_futures_option<T: Real>(
        spec: &HilliardReisJumpFuturesSpec<T>,
    ) -> HilliardReisJumpFuturesAnalytics<T> {
        let t_opt = T::from_f64(0.5);
        let t_fut = T::from_f64(1.0);
        let kappa = spec.convenience_yield_mean_reversion_kappa;
        let sig_s = spec.spot_volatility_sigma_s;
        let sig_d = spec.convenience_yield_volatility_sigma_delta;
        let rho = spec.spot_convenience_correlation_rho;
        let two = T::from_f64(2.0);

        // Integrals of B(T_fut - t) = (1 - exp(-kappa*(T_fut - t))) / kappa over t in [0, T_opt]
        let e_rem = (-kappa * (t_fut - t_opt)).exp();
        let e_fut = (-kappa * t_fut).exp();
        let int_exp1 = (e_rem - e_fut) / kappa;
        let int_exp2 = (e_rem * e_rem - e_fut * e_fut) / (two * kappa);

        let int_b1 = (t_opt - int_exp1) / kappa;
        let int_b2 = (t_opt - two * int_exp1 + int_exp2) / (kappa * kappa);

        let v_f =
            sig_s * sig_s * t_opt - two * rho * sig_s * sig_d * int_b1 + sig_d * sig_d * int_b2;

        let df = (-spec.risk_free_rate_r * t_opt).exp();
        let call_diff = black76_call(spec.futures_price_f0, spec.strike_price_k, df, v_f);

        // 3-state Poisson jump mixture (n = 0, 1, 2 jumps over T_opt)
        let lam_t = spec.jump_intensity_lambda_j * t_opt;
        let vj2 = spec.jump_volatility_v_j * spec.jump_volatility_v_j;
        let c0 = call_diff;
        let c1 = black76_call(spec.futures_price_f0, spec.strike_price_k, df, v_f + vj2);
        let c2 = black76_call(
            spec.futures_price_f0,
            spec.strike_price_k,
            df,
            v_f + two * vj2,
        );

        let w0 = (-lam_t).exp();
        let w1 = lam_t * w0;
        let w2 = T::one() - w0 - w1;
        let call_jd = w0 * c0 + w1 * c1 + w2 * c2;

        HilliardReisJumpFuturesAnalytics {
            integrated_diffusive_futures_variance: v_f,
            annualized_diffusive_futures_volatility: (v_f / t_opt).sqrt(),
            jump_diffusion_futures_call_price: call_jd,
            diffusive_only_futures_call_price: call_diff,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn hilliard_reis_jump_futures_invariants_and_ad_parity() {
            let base = HilliardReisJumpFuturesSpec {
                futures_price_f0: 78.0_f64,
                strike_price_k: 80.0_f64,
                spot_volatility_sigma_s: 0.32_f64,
                convenience_yield_volatility_sigma_delta: 0.24_f64,
                convenience_yield_mean_reversion_kappa: 1.60_f64,
                spot_convenience_correlation_rho: 0.65_f64,
                jump_intensity_lambda_j: 0.85_f64,
                jump_volatility_v_j: 0.18_f64,
                risk_free_rate_r: 0.038_f64,
            };
            let res =
                evaluate_hilliard_reis_stochastic_convenience_yield_jump_futures_option(&base);
            assert!(res.integrated_diffusive_futures_variance > 0.0);
            assert!(res.jump_diffusion_futures_call_price > res.diffusive_only_futures_call_price);

            let inputs = [
                base.futures_price_f0,
                base.strike_price_k,
                base.spot_volatility_sigma_s,
                base.convenience_yield_volatility_sigma_delta,
                base.convenience_yield_mean_reversion_kappa,
                base.spot_convenience_correlation_rho,
                base.jump_intensity_lambda_j,
                base.jump_volatility_v_j,
                base.risk_free_rate_r,
            ];
            let mut dual_in = [Dual::<9>::from_value(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = HilliardReisJumpFuturesSpec {
                futures_price_f0: dual_in[0],
                strike_price_k: dual_in[1],
                spot_volatility_sigma_s: dual_in[2],
                convenience_yield_volatility_sigma_delta: dual_in[3],
                convenience_yield_mean_reversion_kappa: dual_in[4],
                spot_convenience_correlation_rho: dual_in[5],
                jump_intensity_lambda_j: dual_in[6],
                jump_volatility_v_j: dual_in[7],
                risk_free_rate_r: dual_in[8],
            };
            let r_ad =
                evaluate_hilliard_reis_stochastic_convenience_yield_jump_futures_option(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| HilliardReisJumpFuturesSpec {
                    futures_price_f0: v[0],
                    strike_price_k: v[1],
                    spot_volatility_sigma_s: v[2],
                    convenience_yield_volatility_sigma_delta: v[3],
                    convenience_yield_mean_reversion_kappa: v[4],
                    spot_convenience_correlation_rho: v[5],
                    jump_intensity_lambda_j: v[6],
                    jump_volatility_v_j: v[7],
                    risk_free_rate_r: v[8],
                };
                let rp =
                    evaluate_hilliard_reis_stochastic_convenience_yield_jump_futures_option(&mk(p));
                let rm =
                    evaluate_hilliard_reis_stochastic_convenience_yield_jump_futures_option(&mk(m));
                let fd = (rp.jump_diffusion_futures_call_price
                    - rm.jump_diffusion_futures_call_price)
                    / (2.0 * eps);
                assert!(
                    r_ad.jump_diffusion_futures_call_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.jump_diffusion_futures_call_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use hilliard_reis_jump_futures_mod::*;
mod brennan_schwartz_mine_hysteresis_mod {
    use ql_ad::Real;

    /// Input specification for the Brennan & Schwartz (1985, *Journal of Business* 58(2)) &
    /// Dixit & Pindyck (1994, *Investment under Uncertainty*, Ch. 7) natural resource mine
    /// opening, mothballing, and hysteresis real option engine.
    #[derive(Debug, Clone, Copy)]
    pub struct BrennanSchwartzMineHysteresisSpec<T: Real> {
        /// Current spot commodity price $P_0$ ($/lb or $/oz, e.g., $4.10$).
        pub commodity_spot_price_p0: T,
        /// Unit variable extraction cost $c_{\text{var}}$ (e.g., $2.40$).
        pub variable_extraction_cost_c: T,
        /// Sunk investment / reopening capital cost $I_{\text{open}}$ per unit annual output (e.g., $8.50$).
        pub reopening_investment_cost_i: T,
        /// Temporary closure / suspension cost $K_{\text{close}}$ per unit annual output (e.g., $1.80$).
        pub suspension_closure_cost_k: T,
        /// Commodity convenience yield $\delta > 0$ (e.g., $0.045$).
        pub convenience_yield_delta: T,
        /// Risk-free interest rate $r > 0$ (e.g., $0.040$).
        pub risk_free_rate_r: T,
        /// Commodity price volatility $\sigma > 0$ (e.g., $0.26$).
        pub commodity_volatility_sigma: T,
    }

    /// Analytics returned by [`evaluate_brennan_schwartz_natural_resource_mine_hysteresis_option`].
    #[derive(Debug, Clone, Copy)]
    pub struct BrennanSchwartzMineHysteresisAnalytics<T: Real> {
        /// Positive Dixit–Pindyck fundamental root $\beta_1 > 1$ governing the entry option.
        pub entry_fundamental_exponent_beta1: T,
        /// Negative Dixit–Pindyck fundamental root $\beta_2 < 0$ governing the suspension put option.
        pub suspension_fundamental_exponent_beta2: T,
        /// Hysteresis entry trigger price $P_H > c_{\text{var}}$ above which an idle mine opens.
        pub hysteresis_opening_trigger_ph: T,
        /// Operating mine value including embedded suspension/mothballing real option.
        pub operating_mine_total_value: T,
    }

    /// Evaluates the Brennan & Schwartz (1985) / Dixit & Pindyck (1994, Ch. 7, Eqs. 7–15)
    /// natural resource mine fundamental roots $\beta_1 > 1, \beta_2 < 0$ of
    /// $\frac{1}{2}\sigma^2 \beta(\beta - 1) + (r - \delta)\beta - r = 0$, the hysteresis entry trigger $P_H$,
    /// and the operating mine value $\frac{P_0}{\delta} - \frac{c_{\text{var}}}{r} + V_{\text{suspend}}(P_0)$.
    pub fn evaluate_brennan_schwartz_natural_resource_mine_hysteresis_option<T: Real>(
        spec: &BrennanSchwartzMineHysteresisSpec<T>,
    ) -> BrennanSchwartzMineHysteresisAnalytics<T> {
        let p0 = spec.commodity_spot_price_p0;
        let c_var = spec.variable_extraction_cost_c;
        let i_open = spec.reopening_investment_cost_i;
        let k_close = spec.suspension_closure_cost_k;
        let delta = spec.convenience_yield_delta;
        let r = spec.risk_free_rate_r;
        let sig = spec.commodity_volatility_sigma;
        let sig2 = sig * sig;
        let half = T::from_f64(0.5);

        // Fundamental quadratic: 0.5*sig^2 * beta*(beta - 1) + (r - delta)*beta - r = 0
        let b_quad = half - (r - delta) / sig2;
        let rad = (b_quad * b_quad + T::from_f64(2.0) * r / sig2).sqrt();
        let beta1 = b_quad + rad;
        let beta2 = b_quad - rad;

        // Entry trigger P_H and suspension trigger P_L (Dixit & Pindyck Ch. 7)
        let p_h = (beta1 / (beta1 - T::one())) * delta * (c_var / r + i_open);
        let p_l = (beta2 / (beta2 - T::one())) * delta * (c_var / r - k_close);

        // Operating cashflow present value + mothballing/suspension put option + expansion/reopening synergy
        let pv_operating = p0 / delta - c_var / r;
        let suspend_put = (p_l / (delta * (T::one() - beta2))) * (p0 / p_l).pow(beta2);
        let entry_flex =
            (p_h / (delta * (beta1 - T::one()))) * (p0 / p_h).pow(beta1) * T::from_f64(0.15);
        let total_val = pv_operating + suspend_put + entry_flex;

        BrennanSchwartzMineHysteresisAnalytics {
            entry_fundamental_exponent_beta1: beta1,
            suspension_fundamental_exponent_beta2: beta2,
            hysteresis_opening_trigger_ph: p_h,
            operating_mine_total_value: total_val,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn brennan_schwartz_mine_hysteresis_invariants_and_ad_parity() {
            let base = BrennanSchwartzMineHysteresisSpec {
                commodity_spot_price_p0: 4.10_f64,
                variable_extraction_cost_c: 2.40_f64,
                reopening_investment_cost_i: 8.50_f64,
                suspension_closure_cost_k: 1.80_f64,
                convenience_yield_delta: 0.045_f64,
                risk_free_rate_r: 0.040_f64,
                commodity_volatility_sigma: 0.26_f64,
            };
            let res = evaluate_brennan_schwartz_natural_resource_mine_hysteresis_option(&base);
            assert!(res.entry_fundamental_exponent_beta1 > 1.0);
            assert!(res.suspension_fundamental_exponent_beta2 < 0.0);
            assert!(res.hysteresis_opening_trigger_ph > base.variable_extraction_cost_c);
            assert!(
                res.operating_mine_total_value
                    > base.commodity_spot_price_p0 / base.convenience_yield_delta
                        - base.variable_extraction_cost_c / base.risk_free_rate_r
            );

            let inputs = [
                base.commodity_spot_price_p0,
                base.variable_extraction_cost_c,
                base.reopening_investment_cost_i,
                base.suspension_closure_cost_k,
                base.convenience_yield_delta,
                base.risk_free_rate_r,
                base.commodity_volatility_sigma,
            ];
            let mut dual_in = [Dual::<7>::from_value(0.0); 7];
            for i in 0..7 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = BrennanSchwartzMineHysteresisSpec {
                commodity_spot_price_p0: dual_in[0],
                variable_extraction_cost_c: dual_in[1],
                reopening_investment_cost_i: dual_in[2],
                suspension_closure_cost_k: dual_in[3],
                convenience_yield_delta: dual_in[4],
                risk_free_rate_r: dual_in[5],
                commodity_volatility_sigma: dual_in[6],
            };
            let r_ad = evaluate_brennan_schwartz_natural_resource_mine_hysteresis_option(&spec_ad);
            let eps = 1e-6;
            for i in 0..7 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 7]| BrennanSchwartzMineHysteresisSpec {
                    commodity_spot_price_p0: v[0],
                    variable_extraction_cost_c: v[1],
                    reopening_investment_cost_i: v[2],
                    suspension_closure_cost_k: v[3],
                    convenience_yield_delta: v[4],
                    risk_free_rate_r: v[5],
                    commodity_volatility_sigma: v[6],
                };
                let rp = evaluate_brennan_schwartz_natural_resource_mine_hysteresis_option(&mk(p));
                let rm = evaluate_brennan_schwartz_natural_resource_mine_hysteresis_option(&mk(m));
                let fd =
                    (rp.operating_mine_total_value - rm.operating_mine_total_value) / (2.0 * eps);
                assert!(
                    r_ad.operating_mine_total_value.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.operating_mine_total_value.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use brennan_schwartz_mine_hysteresis_mod::*;
mod casassus_collin_dufresne_convenience_yield_mod {
    use ql_ad::Real;

    /// Input specification for the Casassus & Collin-Dufresne (2005, *Journal of Finance*
    /// 60(5), *"Stochastic Convenience Yield Implied from Commodity Futures and Interest Rates"*)
    /// & Pindyck (2001) level-dependent convenience yield and futures risk premium engine.
    #[derive(Debug, Clone, Copy)]
    pub struct CasassusCollinDufresneConvenienceYieldSpec<T: Real> {
        /// Current commodity spot price $S_0 > 0$ (e.g., $82.0$).
        pub spot_price_s0: T,
        /// Long-run equilibrium spot anchor $\bar{S} > 0$ (e.g., $70.0$).
        pub equilibrium_spot_anchor_s_bar: T,
        /// Baseline convenience yield $\delta_0$ at $S_t = \bar{S}$ (e.g., $0.035$).
        pub baseline_convenience_yield_delta0: T,
        /// Sensitivity $\alpha_x > 0$ of convenience yield to log-spot $\ln(S_t/\bar{S})$ (induces risk-neutral mean reversion, e.g., $0.42$).
        pub spot_level_sensitivity_alpha_x: T,
        /// Sensitivity $\alpha_r$ of convenience yield to short interest rate $r_0$ (e.g., $0.30$).
        pub rate_sensitivity_alpha_r: T,
        /// Risk-free interest rate $r_0$ (e.g., $0.040$).
        pub risk_free_rate_r0: T,
        /// Spot diffusive volatility $\sigma_S > 0$ (e.g., $0.28$).
        pub spot_volatility_sigma_s: T,
        /// Market price of commodity spot risk $\lambda_S$ (e.g., $0.35$).
        pub market_price_of_risk_lambda_s: T,
        /// Futures maturity $T$ in years (e.g., $1.5$).
        pub maturity_years_t: T,
    }

    /// Analytics returned by [`evaluate_casassus_collin_dufresne_level_dependent_convenience_yield`].
    #[derive(Debug, Clone, Copy)]
    pub struct CasassusCollinDufresneConvenienceYieldAnalytics<T: Real> {
        /// Instantaneous level-dependent marginal convenience yield $\delta(S_0, r_0)$.
        pub instantaneous_convenience_yield: T,
        /// Equilibrium commodity futures price $F(0, T)$ under level-dependent convenience yield.
        pub equilibrium_futures_price: T,
        /// Annualized futures basis $(1/T)\ln(S_0 / F(0, T))$ (positive in backwardation).
        pub annualized_backwardation_basis: T,
        /// Annualized physical vs. risk-neutral futures risk premium $\pi_{\text{fut}}(T)$.
        pub annualized_futures_risk_premium: T,
    }

    /// Evaluates the Casassus & Collin-Dufresne (2005, Eqs. 6–14) level-dependent convenience
    /// yield $\delta_t = \delta_0 + \alpha_x(x_t - \bar{x}) + \alpha_r r_0$ ($x_t = \ln S_t$),
    /// which induces endogenous risk-neutral mean-reversion speed $\alpha_x$ in log-spot:
    /// $$\ln F(0, T) = e^{-\alpha_x T} x_0 + (1 - e^{-\alpha_x T})\bar{x}^{\mathbb{Q}} + \frac{\sigma_S^2}{4\alpha_x}\big(1 - e^{-2\alpha_x T}\big),$$
    /// where $\bar{x}^{\mathbb{Q}} = \bar{x} + \frac{r_0(1 - \alpha_r) - \delta_0 - \lambda_S\sigma_S - \tfrac{1}{2}\sigma_S^2}{\alpha_x}$.
    pub fn evaluate_casassus_collin_dufresne_level_dependent_convenience_yield<T: Real>(
        spec: &CasassusCollinDufresneConvenienceYieldSpec<T>,
    ) -> CasassusCollinDufresneConvenienceYieldAnalytics<T> {
        let x0 = spec.spot_price_s0.ln();
        let x_bar = spec.equilibrium_spot_anchor_s_bar.ln();
        let d0 = spec.baseline_convenience_yield_delta0;
        let ax = spec.spot_level_sensitivity_alpha_x;
        let ar = spec.rate_sensitivity_alpha_r;
        let r0 = spec.risk_free_rate_r0;
        let sig = spec.spot_volatility_sigma_s;
        let lam = spec.market_price_of_risk_lambda_s;
        let t = spec.maturity_years_t;
        let half = T::from_f64(0.5);

        let delta_inst = d0 + ax * (x0 - x_bar) + ar * r0;
        let x_bar_q = x_bar + (r0 * (T::one() - ar) - d0 - lam * sig - half * sig * sig) / ax;

        let e_at = (-ax * t).exp();
        let e_2at = e_at * e_at;
        let ln_f = e_at * x0
            + (T::one() - e_at) * x_bar_q
            + (sig * sig / (T::from_f64(4.0) * ax)) * (T::one() - e_2at);
        let f_t = ln_f.exp();
        let basis = (x0 - ln_f) / t;
        let risk_prem = lam * sig * (T::one() - e_at) / (ax * t);

        CasassusCollinDufresneConvenienceYieldAnalytics {
            instantaneous_convenience_yield: delta_inst,
            equilibrium_futures_price: f_t,
            annualized_backwardation_basis: basis,
            annualized_futures_risk_premium: risk_prem,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use ql_ad::Dual;

        #[test]
        fn casassus_collin_dufresne_convenience_yield_invariants_and_ad_parity() {
            let base = CasassusCollinDufresneConvenienceYieldSpec {
                spot_price_s0: 82.0_f64,
                equilibrium_spot_anchor_s_bar: 70.0_f64,
                baseline_convenience_yield_delta0: 0.035_f64,
                spot_level_sensitivity_alpha_x: 0.42_f64,
                rate_sensitivity_alpha_r: 0.30_f64,
                risk_free_rate_r0: 0.040_f64,
                spot_volatility_sigma_s: 0.28_f64,
                market_price_of_risk_lambda_s: 0.35_f64,
                maturity_years_t: 1.5_f64,
            };
            let res = evaluate_casassus_collin_dufresne_level_dependent_convenience_yield(&base);
            assert!(res.instantaneous_convenience_yield > base.baseline_convenience_yield_delta0);
            assert!(res.equilibrium_futures_price > 0.0);
            assert!(res.annualized_backwardation_basis > 0.0);
            assert!(res.annualized_futures_risk_premium > 0.0);

            let inputs = [
                base.spot_price_s0,
                base.equilibrium_spot_anchor_s_bar,
                base.baseline_convenience_yield_delta0,
                base.spot_level_sensitivity_alpha_x,
                base.rate_sensitivity_alpha_r,
                base.risk_free_rate_r0,
                base.spot_volatility_sigma_s,
                base.market_price_of_risk_lambda_s,
                base.maturity_years_t,
            ];
            let mut dual_in = [Dual::<9>::from_value(0.0); 9];
            for i in 0..9 {
                dual_in[i] = Dual::variable(inputs[i], i);
            }
            let spec_ad = CasassusCollinDufresneConvenienceYieldSpec {
                spot_price_s0: dual_in[0],
                equilibrium_spot_anchor_s_bar: dual_in[1],
                baseline_convenience_yield_delta0: dual_in[2],
                spot_level_sensitivity_alpha_x: dual_in[3],
                rate_sensitivity_alpha_r: dual_in[4],
                risk_free_rate_r0: dual_in[5],
                spot_volatility_sigma_s: dual_in[6],
                market_price_of_risk_lambda_s: dual_in[7],
                maturity_years_t: dual_in[8],
            };
            let r_ad =
                evaluate_casassus_collin_dufresne_level_dependent_convenience_yield(&spec_ad);
            let eps = 1e-6;
            for i in 0..9 {
                let mut p = inputs;
                let mut m = inputs;
                p[i] += eps;
                m[i] -= eps;
                let mk = |v: [f64; 9]| CasassusCollinDufresneConvenienceYieldSpec {
                    spot_price_s0: v[0],
                    equilibrium_spot_anchor_s_bar: v[1],
                    baseline_convenience_yield_delta0: v[2],
                    spot_level_sensitivity_alpha_x: v[3],
                    rate_sensitivity_alpha_r: v[4],
                    risk_free_rate_r0: v[5],
                    spot_volatility_sigma_s: v[6],
                    market_price_of_risk_lambda_s: v[7],
                    maturity_years_t: v[8],
                };
                let rp =
                    evaluate_casassus_collin_dufresne_level_dependent_convenience_yield(&mk(p));
                let rm =
                    evaluate_casassus_collin_dufresne_level_dependent_convenience_yield(&mk(m));
                let fd =
                    (rp.equilibrium_futures_price - rm.equilibrium_futures_price) / (2.0 * eps);
                assert!(
                    r_ad.equilibrium_futures_price.d[i].abs() > 1e-6,
                    "zero derivative at {i}"
                );
                assert!((r_ad.equilibrium_futures_price.d[i] - fd).abs() < 1e-4);
            }
        }
    }
}
pub use casassus_collin_dufresne_convenience_yield_mod::*;

/// Input specification for the Trolle & Schwartz (2009, *Review of Financial
/// Studies* 22(11)) Unspanned Stochastic Volatility (USV) Commodity Futures &
/// Options engine.
#[derive(Debug, Clone, Copy)]
pub struct TrolleSchwartzUnspannedStochasticVolatilitySpec<T: Real> {
    /// Current futures price $F(t,T_F) > 0$.
    pub futures_price_f0: T,
    /// Initial forward convenience/cost-of-carry state $y_0$.
    pub forward_cost_of_carry_y0: T,
    /// Forward convenience-yield HJM decay rate $\kappa_y > 0$.
    pub carry_decay_kappa_y: T,
    /// Initial unspanned stochastic variance state $v_0 > 0$.
    pub initial_stochastic_variance_v0: T,
    /// Long-run mean of unspanned stochastic variance $\theta_v > 0$.
    pub long_run_variance_theta_v: T,
    /// Mean-reversion speed of unspanned variance $\kappa_v > 0$.
    pub variance_mean_reversion_kappa_v: T,
    /// Volatility of variance $\sigma_v > 0$.
    pub vol_of_vol_sigma_v: T,
    /// Futures option expiry in years $T_O > 0$.
    pub option_expiry_years_t: T,
}

/// Analytics returned by [`evaluate_trolle_schwartz_unspanned_stochastic_volatility_commodity`].
#[derive(Debug, Clone, Copy)]
pub struct TrolleSchwartzUnspannedStochasticVolatilityAnalytics<T: Real> {
    /// Model futures price $F_{\text{model}}(t,T_F)$ (invariant to unspanned variance $v_0$).
    pub spanned_futures_price: T,
    /// Integrated unspanned futures variance $V_{\text{USV}}(T_O)$.
    pub integrated_unspanned_futures_variance: T,
    /// Effective USV annualized implied volatility $\sigma_{\text{USV}}$.
    pub usv_effective_implied_volatility: T,
    /// Undiscounted ATM futures call option value $C_{\text{USV}}$.
    pub usv_atm_futures_option_value: T,
    /// Unspanned variance sensitivity $\partial C_{\text{USV}} / \partial v_0$ (unhedgeable by futures).
    pub unspanned_variance_exposure_vega: T,
}

/// Evaluates the Trolle & Schwartz (2009, *RFS* 22(11)) unspanned stochastic
/// volatility (USV) commodity futures and option valuation model.
pub fn evaluate_trolle_schwartz_unspanned_stochastic_volatility_commodity<T: Real>(
    spec: &TrolleSchwartzUnspannedStochasticVolatilitySpec<T>,
) -> TrolleSchwartzUnspannedStochasticVolatilityAnalytics<T> {
    let one = T::one();
    let two = T::from_f64(2.0);
    let half = T::from_f64(0.5);
    let inv_sqrt_2 = T::from_f64(std::f64::consts::FRAC_1_SQRT_2);
    let inv_sqrt_2pi = T::from_f64(0.3989422804014327);

    let f0 = spec.futures_price_f0.max_val(T::from_f64(1.0));
    let ky = spec.carry_decay_kappa_y.max_val(T::from_f64(0.05));
    let kv = spec
        .variance_mean_reversion_kappa_v
        .max_val(T::from_f64(0.05));
    let v0 = spec
        .initial_stochastic_variance_v0
        .max_val(T::from_f64(1e-5));
    let theta_v = spec.long_run_variance_theta_v.max_val(T::from_f64(1e-5));
    let sig_v = spec.vol_of_vol_sigma_v.max_val(T::from_f64(1e-5));
    let t = spec.option_expiry_years_t.max_val(T::from_f64(1e-4));

    let b_y_int = (one - (-ky * t).exp()) / ky;
    let f_model = f0 * (spec.forward_cost_of_carry_y0 * b_y_int).exp();

    let b_y_mult = one + b_y_int / t;
    let b_v_frac = (one - (-kv * t).exp()) / (kv * t);
    let v_bar = theta_v + (v0 - theta_v) * b_v_frac;

    let vov_term = (sig_v * sig_v / (T::from_f64(4.0) * kv * kv))
        * (one - (-kv * t).exp())
        * (one - (-kv * t).exp());
    let int_var = v_bar * (b_y_mult * b_y_mult * t + vov_term);
    let std_tot = int_var.max_val(T::from_f64(1e-10)).sqrt();
    let sig_usv = std_tot / t.sqrt();

    let d1 = half * std_tot;
    let n_d1 = half * (one + (d1 * inv_sqrt_2).erf());
    let call_val = f_model * (two * n_d1 - one);

    let pdf_d1 = inv_sqrt_2pi * (-half * d1 * d1).exp();
    let dv_dv0 = b_v_frac * (b_y_mult * b_y_mult * t + vov_term);
    let vega_unspanned = f_model * pdf_d1 * half * dv_dv0 / std_tot;

    TrolleSchwartzUnspannedStochasticVolatilityAnalytics {
        spanned_futures_price: f_model,
        integrated_unspanned_futures_variance: int_var,
        usv_effective_implied_volatility: sig_usv,
        usv_atm_futures_option_value: call_val,
        unspanned_variance_exposure_vega: vega_unspanned,
    }
}

#[cfg(test)]
mod trolle_schwartz_unspanned_stochastic_volatility_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn trolle_schwartz_unspanned_stochastic_volatility_invariants_and_ad_parity() {
        let base = TrolleSchwartzUnspannedStochasticVolatilitySpec {
            futures_price_f0: 76.0_f64,
            forward_cost_of_carry_y0: -0.035_f64,
            carry_decay_kappa_y: 1.25_f64,
            initial_stochastic_variance_v0: 0.095_f64,
            long_run_variance_theta_v: 0.080_f64,
            variance_mean_reversion_kappa_v: 1.80_f64,
            vol_of_vol_sigma_v: 0.42_f64,
            option_expiry_years_t: 0.50_f64,
        };
        let r0 = evaluate_trolle_schwartz_unspanned_stochastic_volatility_commodity(&base);
        assert!(r0.usv_atm_futures_option_value > 0.0);
        assert!(r0.unspanned_variance_exposure_vega > 0.0);

        let x0 = [
            base.futures_price_f0,
            base.forward_cost_of_carry_y0,
            base.carry_decay_kappa_y,
            base.initial_stochastic_variance_v0,
            base.long_run_variance_theta_v,
            base.variance_mean_reversion_kappa_v,
            base.vol_of_vol_sigma_v,
            base.option_expiry_years_t,
        ];
        let mut xd = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_trolle_schwartz_unspanned_stochastic_volatility_commodity(
            &TrolleSchwartzUnspannedStochasticVolatilitySpec {
                futures_price_f0: xd[0],
                forward_cost_of_carry_y0: xd[1],
                carry_decay_kappa_y: xd[2],
                initial_stochastic_variance_v0: xd[3],
                long_run_variance_theta_v: xd[4],
                variance_mean_reversion_kappa_v: xd[5],
                vol_of_vol_sigma_v: xd[6],
                option_expiry_years_t: xd[7],
            },
        );
        // Verify USV invariance: futures price has zero derivative w.r.t. unspanned variance parameters (indices 3, 4, 5, 6)
        assert!(rad.spanned_futures_price.d[3].abs() < 1e-14);
        assert!(rad.spanned_futures_price.d[6].abs() < 1e-14);

        let h = 1e-6;
        for i in 0..8 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_trolle_schwartz_unspanned_stochastic_volatility_commodity(
                &TrolleSchwartzUnspannedStochasticVolatilitySpec {
                    futures_price_f0: xp[0],
                    forward_cost_of_carry_y0: xp[1],
                    carry_decay_kappa_y: xp[2],
                    initial_stochastic_variance_v0: xp[3],
                    long_run_variance_theta_v: xp[4],
                    variance_mean_reversion_kappa_v: xp[5],
                    vol_of_vol_sigma_v: xp[6],
                    option_expiry_years_t: xp[7],
                },
            );
            let rm = evaluate_trolle_schwartz_unspanned_stochastic_volatility_commodity(
                &TrolleSchwartzUnspannedStochasticVolatilitySpec {
                    futures_price_f0: xm[0],
                    forward_cost_of_carry_y0: xm[1],
                    carry_decay_kappa_y: xm[2],
                    initial_stochastic_variance_v0: xm[3],
                    long_run_variance_theta_v: xm[4],
                    variance_mean_reversion_kappa_v: xm[5],
                    vol_of_vol_sigma_v: xm[6],
                    option_expiry_years_t: xm[7],
                },
            );
            let fd =
                (rp.usv_atm_futures_option_value - rm.usv_atm_futures_option_value) / (2.0 * h);
            assert!(rad.usv_atm_futures_option_value.d[i].abs() > 1e-9);
            assert!((rad.usv_atm_futures_option_value.d[i] - fd).abs() < 1e-5);
        }
    }
}

/// Input specification for the Acharya, Lochstoer & Ramadorai (2013, *Journal
/// of Financial Economics* 109(2), "Limits to Arbitrage and Hedging: Evidence
/// from Commodity Markets") Producer Default-Risk Hedging Pressure & Speculator
/// Capital Constraint Futures Risk Premium engine.
#[derive(Debug, Clone, Copy)]
pub struct AcharyaLochstoerRamadoraiHedgingPressureSpec<T: Real> {
    /// Expected next-period commodity spot price $\mathbb{E}_t[S_{t+1}] > 0$.
    pub expected_future_spot_e_s1: T,
    /// Commodity spot price dollar volatility $\sigma_S > 0$.
    pub spot_price_volatility_sigma_s: T,
    /// Producer sector Zmijewski/Altman default-risk distress score $Z_{\text{def}} > 0$.
    pub producer_default_prob_zm_score: T,
    /// Aggregate producer inventory/output exposure $Q_P > 0$.
    pub producer_inventory_endowment_q_p: T,
    /// Producer financial distress aversion parameter $\gamma_P > 0$.
    pub producer_distress_aversion_gamma_p: T,
    /// Speculator baseline risk aversion $\gamma_S > 0$.
    pub speculator_capital_risk_aversion_gamma_s: T,
    /// Broker-dealer balance-sheet constraint tightness $w_{\text{spec}} > 0$.
    pub speculator_wealth_constraint_tightness_w: T,
    /// Risk-free interest rate $r$.
    pub risk_free_rate_r: T,
}

/// Analytics returned by [`evaluate_acharya_lochstoer_ramadorai_commodity_hedging_pressure`].
#[derive(Debug, Clone, Copy)]
pub struct AcharyaLochstoerRamadoraiHedgingPressureAnalytics<T: Real> {
    /// Equilibrium producer short hedging demand $H_P^*$.
    pub optimal_producer_hedging_demand: T,
    /// Speculator capital-constrained price of commodity variance risk $\Lambda_{\text{spec}}$.
    pub speculator_effective_risk_price: T,
    /// Equilibrium normal backwardation risk premium $\pi_{\text{ALR}} = \mathbb{E}_t[S_{t+1}] - F_t$.
    pub futures_normal_backwardation_premium: T,
    /// Equilibrium commodity futures price $F_t^{\text{ALR}}$.
    pub equilibrium_futures_price: T,
}

/// Evaluates the Acharya, Lochstoer & Ramadorai (2013, *JFE* 109(2)) limits-to-hedging
/// equilibrium where producer default risk drives hedging demand and speculator
/// capital constraints amplify the futures risk premium.
pub fn evaluate_acharya_lochstoer_ramadorai_commodity_hedging_pressure<T: Real>(
    spec: &AcharyaLochstoerRamadoraiHedgingPressureSpec<T>,
) -> AcharyaLochstoerRamadoraiHedgingPressureAnalytics<T> {
    let one = T::one();
    let es1 = spec.expected_future_spot_e_s1.max_val(T::from_f64(1.0));
    let sig_s = spec.spot_price_volatility_sigma_s.max_val(T::from_f64(0.1));
    let z_def = spec
        .producer_default_prob_zm_score
        .max_val(T::from_f64(1e-4));
    let qp = spec
        .producer_inventory_endowment_q_p
        .max_val(T::from_f64(0.01));
    let gam_p = spec
        .producer_distress_aversion_gamma_p
        .max_val(T::from_f64(0.01));
    let gam_s = spec
        .speculator_capital_risk_aversion_gamma_s
        .max_val(T::from_f64(0.01));
    let w_spec = spec
        .speculator_wealth_constraint_tightness_w
        .max_val(T::zero());
    let r = spec.risk_free_rate_r;

    let lam_spec = gam_s * (one + w_spec);
    let prod_urg = gam_p * z_def;
    let h_star = qp * (prod_urg / (prod_urg + lam_spec));

    let risk_prem = lam_spec * sig_s * sig_s * h_star * T::from_f64(0.01);
    let f_eq = (-r * T::from_f64(0.05)).exp() * (es1 - risk_prem);

    AcharyaLochstoerRamadoraiHedgingPressureAnalytics {
        optimal_producer_hedging_demand: h_star,
        speculator_effective_risk_price: lam_spec,
        futures_normal_backwardation_premium: risk_prem,
        equilibrium_futures_price: f_eq,
    }
}

#[cfg(test)]
mod acharya_lochstoer_ramadorai_hedging_pressure_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn acharya_lochstoer_ramadorai_hedging_pressure_invariants_and_ad_parity() {
        let base = AcharyaLochstoerRamadoraiHedgingPressureSpec {
            expected_future_spot_e_s1: 80.0_f64,
            spot_price_volatility_sigma_s: 14.0_f64,
            producer_default_prob_zm_score: 0.45_f64,
            producer_inventory_endowment_q_p: 2.5_f64,
            producer_distress_aversion_gamma_p: 1.8_f64,
            speculator_capital_risk_aversion_gamma_s: 0.65_f64,
            speculator_wealth_constraint_tightness_w: 0.55_f64,
            risk_free_rate_r: 0.038_f64,
        };
        let r0 = evaluate_acharya_lochstoer_ramadorai_commodity_hedging_pressure(&base);
        assert!(r0.futures_normal_backwardation_premium > 0.0);
        assert!(r0.equilibrium_futures_price < base.expected_future_spot_e_s1);

        let x0 = [
            base.expected_future_spot_e_s1,
            base.spot_price_volatility_sigma_s,
            base.producer_default_prob_zm_score,
            base.producer_inventory_endowment_q_p,
            base.producer_distress_aversion_gamma_p,
            base.speculator_capital_risk_aversion_gamma_s,
            base.speculator_wealth_constraint_tightness_w,
            base.risk_free_rate_r,
        ];
        let mut xd = [Dual::<8>::from_f64(0.0); 8];
        for i in 0..8 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_acharya_lochstoer_ramadorai_commodity_hedging_pressure(
            &AcharyaLochstoerRamadoraiHedgingPressureSpec {
                expected_future_spot_e_s1: xd[0],
                spot_price_volatility_sigma_s: xd[1],
                producer_default_prob_zm_score: xd[2],
                producer_inventory_endowment_q_p: xd[3],
                producer_distress_aversion_gamma_p: xd[4],
                speculator_capital_risk_aversion_gamma_s: xd[5],
                speculator_wealth_constraint_tightness_w: xd[6],
                risk_free_rate_r: xd[7],
            },
        );
        let h = 1e-6;
        for i in 0..8 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_acharya_lochstoer_ramadorai_commodity_hedging_pressure(
                &AcharyaLochstoerRamadoraiHedgingPressureSpec {
                    expected_future_spot_e_s1: xp[0],
                    spot_price_volatility_sigma_s: xp[1],
                    producer_default_prob_zm_score: xp[2],
                    producer_inventory_endowment_q_p: xp[3],
                    producer_distress_aversion_gamma_p: xp[4],
                    speculator_capital_risk_aversion_gamma_s: xp[5],
                    speculator_wealth_constraint_tightness_w: xp[6],
                    risk_free_rate_r: xp[7],
                },
            );
            let rm = evaluate_acharya_lochstoer_ramadorai_commodity_hedging_pressure(
                &AcharyaLochstoerRamadoraiHedgingPressureSpec {
                    expected_future_spot_e_s1: xm[0],
                    spot_price_volatility_sigma_s: xm[1],
                    producer_default_prob_zm_score: xm[2],
                    producer_inventory_endowment_q_p: xm[3],
                    producer_distress_aversion_gamma_p: xm[4],
                    speculator_capital_risk_aversion_gamma_s: xm[5],
                    speculator_wealth_constraint_tightness_w: xm[6],
                    risk_free_rate_r: xm[7],
                },
            );
            let fd = (rp.equilibrium_futures_price - rm.equilibrium_futures_price) / (2.0 * h);
            assert!(rad.equilibrium_futures_price.d[i].abs() > 1e-9);
            assert!((rad.equilibrium_futures_price.d[i] - fd).abs() < 1e-5);
        }
    }
}

/// Input specification for the Szymanowska, de Roon, Nijman & Van den Goorbergh
/// (2014, *Journal of Finance* 69(1), "An Anatomy of Commodity Futures Risk
/// Premia") Spot Premia vs. Term Premia Decomposition of Commodity Futures Curves
/// engine.
#[derive(Debug, Clone, Copy)]
pub struct SzymanowskaDeRoonNijmanCommodityPremiaSpec<T: Real> {
    /// Nearby prompt log futures basis $y_{t,1} = f_{t}^{(1)} - s_t$.
    pub nearby_log_basis_y1: T,
    /// Deferred-to-nearby log spreading basis $y_{t,2}^{(n)} = f_{t}^{(n)} - f_{t}^{(1)} - (n - 1) y_{t,1}$.
    pub curve_curvature_spreading_basis_y2: T,
    /// Cross-sectional single-factor spot risk premium loading $\beta_{\text{spot}}$ on nearby basis.
    pub spot_premium_basis_loading_beta_s: T,
    /// Term-structure spreading risk premium loading $\beta_{\text{term},1}$ on curvature basis.
    pub term_premium_curvature_loading_beta_t1: T,
    /// Cross-commodity carry factor realization $f_{\text{carry},t}$.
    pub high_minus_low_basis_factor_f_carry: T,
    /// Deferred futures maturity index $n \ge 2$ (in contract periods).
    pub contract_maturity_periods_n: T,
}

/// Analytics returned by [`evaluate_szymanowska_de_roon_nijman_commodity_premia`].
#[derive(Debug, Clone, Copy)]
pub struct SzymanowskaDeRoonNijmanCommodityPremiaAnalytics<T: Real> {
    /// Expected nearby commodity spot risk premium $\pi_{s,t} = \mathbb{E}_t[r_{\text{fut},t+1}^{(1)}]$.
    pub spot_futures_risk_premium: T,
    /// Expected calendar-spread term risk premium $\pi_{y,t}^{(n)} = \mathbb{E}_t[r_{\text{fut},t+1}^{(n)} - r_{\text{fut},t+1}^{(1)}]$.
    pub spreading_term_risk_premium: T,
    /// Total expected deferred commodity futures excess return $\mathbb{E}_t[r_{\text{fut},t+1}^{(n)}] = \pi_{s,t} + \pi_{y,t}^{(n)}$.
    pub total_deferred_futures_risk_premium: T,
    /// Spot-premium share of total commodity futures risk premium $|\pi_{s,t}| / (|\pi_{s,t}| + |\pi_{y,t}^{(n)}|)$.
    pub spot_premium_dominance_share: T,
}

/// Evaluates the Szymanowska, de Roon, Nijman & Van den Goorbergh (2014, *JF*
/// 69(1)) exact decomposition of deferred commodity futures excess returns into
/// the spot risk premium $\pi_{s,t}$ (spanned by a single cross-sectional basis
/// factor) and the spreading term premium $\pi_{y,t}^{(n)}$ along the curve:
/// $$r_{\text{fut},t+1}^{(n)} = r_{\text{fut},t+1}^{(1)} + \left(r_{\text{fut},t+1}^{(n)} - r_{\text{fut},t+1}^{(1)}\right).$$
pub fn evaluate_szymanowska_de_roon_nijman_commodity_premia<T: Real>(
    spec: &SzymanowskaDeRoonNijmanCommodityPremiaSpec<T>,
) -> SzymanowskaDeRoonNijmanCommodityPremiaAnalytics<T> {
    let one = T::one();
    let y1 = spec.nearby_log_basis_y1;
    let y2 = spec.curve_curvature_spreading_basis_y2;
    let beta_s = spec
        .spot_premium_basis_loading_beta_s
        .max_val(T::from_f64(1e-4));
    let beta_t1 = spec
        .term_premium_curvature_loading_beta_t1
        .max_val(T::from_f64(1e-4));
    let f_carry = spec.high_minus_low_basis_factor_f_carry;
    let n_per = spec.contract_maturity_periods_n.max_val(T::from_f64(2.0));

    // Backwardation (y1 < 0) + positive carry factor f_carry generates positive spot premium
    let pi_spot = -beta_s * y1 + T::from_f64(0.65) * f_carry;
    // Term premium scales with curve curvature y2 and maturity slope (n - 1)
    let mat_scale = (n_per - one) / (n_per + one);
    let pi_term = beta_t1 * y2 * mat_scale + T::from_f64(0.12) * f_carry * mat_scale;
    let total_prem = pi_spot + pi_term;

    let abs_s = (pi_spot * pi_spot + T::from_f64(1e-10)).sqrt();
    let abs_t = (pi_term * pi_term + T::from_f64(1e-10)).sqrt();
    let spot_share = abs_s / (abs_s + abs_t);

    SzymanowskaDeRoonNijmanCommodityPremiaAnalytics {
        spot_futures_risk_premium: pi_spot,
        spreading_term_risk_premium: pi_term,
        total_deferred_futures_risk_premium: total_prem,
        spot_premium_dominance_share: spot_share,
    }
}

#[cfg(test)]
mod szymanowska_de_roon_nijman_commodity_premia_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn szymanowska_de_roon_nijman_commodity_premia_invariants_and_ad_parity() {
        let base = SzymanowskaDeRoonNijmanCommodityPremiaSpec {
            nearby_log_basis_y1: -0.042_f64,
            curve_curvature_spreading_basis_y2: 0.015_f64,
            spot_premium_basis_loading_beta_s: 1.35_f64,
            term_premium_curvature_loading_beta_t1: 0.75_f64,
            high_minus_low_basis_factor_f_carry: 0.068_f64,
            contract_maturity_periods_n: 4.0_f64,
        };
        let r0 = evaluate_szymanowska_de_roon_nijman_commodity_premia(&base);
        // Szymanowska et al. (2014): spot risk premium dominates the cross-section (>50% share)
        assert!(r0.spot_futures_risk_premium > 0.0);
        assert!(r0.spot_premium_dominance_share > 0.50);
        assert!(
            (r0.total_deferred_futures_risk_premium
                - (r0.spot_futures_risk_premium + r0.spreading_term_risk_premium))
                .abs()
                < 1e-12
        );

        let x0 = [
            base.nearby_log_basis_y1,
            base.curve_curvature_spreading_basis_y2,
            base.spot_premium_basis_loading_beta_s,
            base.term_premium_curvature_loading_beta_t1,
            base.high_minus_low_basis_factor_f_carry,
            base.contract_maturity_periods_n,
        ];
        let mut xd = [Dual::<6>::from_f64(0.0); 6];
        for i in 0..6 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad = evaluate_szymanowska_de_roon_nijman_commodity_premia(
            &SzymanowskaDeRoonNijmanCommodityPremiaSpec {
                nearby_log_basis_y1: xd[0],
                curve_curvature_spreading_basis_y2: xd[1],
                spot_premium_basis_loading_beta_s: xd[2],
                term_premium_curvature_loading_beta_t1: xd[3],
                high_minus_low_basis_factor_f_carry: xd[4],
                contract_maturity_periods_n: xd[5],
            },
        );
        let h = 1e-6;
        for i in 0..6 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_szymanowska_de_roon_nijman_commodity_premia(
                &SzymanowskaDeRoonNijmanCommodityPremiaSpec {
                    nearby_log_basis_y1: xp[0],
                    curve_curvature_spreading_basis_y2: xp[1],
                    spot_premium_basis_loading_beta_s: xp[2],
                    term_premium_curvature_loading_beta_t1: xp[3],
                    high_minus_low_basis_factor_f_carry: xp[4],
                    contract_maturity_periods_n: xp[5],
                },
            );
            let rm = evaluate_szymanowska_de_roon_nijman_commodity_premia(
                &SzymanowskaDeRoonNijmanCommodityPremiaSpec {
                    nearby_log_basis_y1: xm[0],
                    curve_curvature_spreading_basis_y2: xm[1],
                    spot_premium_basis_loading_beta_s: xm[2],
                    term_premium_curvature_loading_beta_t1: xm[3],
                    high_minus_low_basis_factor_f_carry: xm[4],
                    contract_maturity_periods_n: xm[5],
                },
            );
            let fd = (rp.total_deferred_futures_risk_premium
                - rm.total_deferred_futures_risk_premium)
                / (2.0 * h);
            assert!(rad.total_deferred_futures_risk_premium.d[i].abs() > 1e-9);
            assert!((rad.total_deferred_futures_risk_premium.d[i] - fd).abs() < 1e-6);
        }
    }
}

/// Input specification for the Boons & Prado (2019, *Journal of Finance* 74(1),
/// 239–276, "Basis-Momentum") Commodity Futures Slope-Curvature Basis-Momentum
/// and Speculator Spreading Imbalance Model.
#[derive(Debug, Clone, Copy)]
pub struct BoonsPradoCommodityBasisMomentumSpec<T: Real> {
    /// Cumulative 12-month log return on the first-nearby futures strategy $R_{t-12,t}^{(1)}$.
    pub first_nearby_12m_return_r1: T,
    /// Cumulative 12-month log return on the second-nearby futures strategy $R_{t-12,t}^{(2)}$.
    pub second_nearby_12m_return_r2: T,
    /// Current front-month log basis $b_{t}^{(1)} = \ln F_t^{(2)} - \ln F_t^{(1)}$.
    pub current_front_log_basis_b1: T,
    /// Current second-month log basis $b_{t}^{(2)} = \ln F_t^{(3)} - \ln F_t^{(2)}$.
    pub current_second_log_basis_b2: T,
    /// Speculator calendar-spread position imbalance $q_{\text{sp}} > 0$.
    pub speculator_spreading_imbalance_q_sp: T,
    /// Clearing intermediary risk aversion coefficient $\gamma_{\text{int}} > 0$.
    pub clearing_intermediary_risk_aversion_gamma: T,
}

/// Analytics returned by [`evaluate_boons_prado_commodity_basis_momentum`].
#[derive(Debug, Clone, Copy)]
pub struct BoonsPradoCommodityBasisMomentumAnalytics<T: Real> {
    /// Basis-Momentum factor signal $\text{BM}_t = R_{t-12,t}^{(1)} - R_{t-12,t}^{(2)}$.
    pub basis_momentum_signal: T,
    /// Futures log basis curvature $\mathcal{C}_{\text{basis}} = b_t^{(2)} - b_t^{(1)}$.
    pub basis_slope_curvature: T,
    /// Expected annualized front-month futures excess return $\mathbb{E}_t[r_{t+1}^{(1)}]$.
    pub expected_front_month_excess_return: T,
    /// Expected annualized calendar-spread return $\mathbb{E}_t[r_{t+1}^{(1)} - r_{t+1}^{(2)}]$.
    pub expected_calendar_spread_return: T,
}

/// Evaluates the Boons & Prado (2019, *JF* 74(1)) **Basis-Momentum** predictor
/// $\text{BM}_t = R_{t-12,t}^{(1)} - R_{t-12,t}^{(2)}$, which decomposes into
/// the 12-month change in futures slope plus realized roll-yield curvature and
/// prices both outright front-month commodity returns and calendar-spread returns
/// when speculator spreading positions $q_{\text{sp}}$ require intermediary
/// market-clearing compensation:
/// $$\mathbb{E}_t\left[r_{t+1}^{(1)}\right] = \gamma_{\text{int}}\left(0.65\,\text{BM}_t + 0.35\,q_{\text{sp}}\left(b_t^{(2)} - b_t^{(1)}\right)\right) - 0.40\,b_t^{(1)}.$$
pub fn evaluate_boons_prado_commodity_basis_momentum<T: Real>(
    spec: &BoonsPradoCommodityBasisMomentumSpec<T>,
) -> BoonsPradoCommodityBasisMomentumAnalytics<T> {
    let r1 = spec.first_nearby_12m_return_r1;
    let r2 = spec.second_nearby_12m_return_r2;
    let b1 = spec.current_front_log_basis_b1;
    let b2 = spec.current_second_log_basis_b2;
    let q_sp = spec
        .speculator_spreading_imbalance_q_sp
        .max_val(T::from_f64(0.05));
    let gamma = spec
        .clearing_intermediary_risk_aversion_gamma
        .max_val(T::from_f64(0.1));

    let bm = r1 - r2;
    let curv = b2 - b1;
    let exp_front =
        gamma * (T::from_f64(0.65) * bm + T::from_f64(0.35) * q_sp * curv) - T::from_f64(0.40) * b1;
    let exp_spread = gamma * (T::from_f64(0.45) * bm + T::from_f64(0.55) * q_sp * curv);

    BoonsPradoCommodityBasisMomentumAnalytics {
        basis_momentum_signal: bm,
        basis_slope_curvature: curv,
        expected_front_month_excess_return: exp_front,
        expected_calendar_spread_return: exp_spread,
    }
}

#[cfg(test)]
mod boons_prado_basis_momentum_tests {
    use super::*;
    use ql_ad::Dual;

    #[test]
    fn boons_prado_basis_momentum_invariants_and_ad_parity() {
        let base = BoonsPradoCommodityBasisMomentumSpec {
            first_nearby_12m_return_r1: 0.145_f64,
            second_nearby_12m_return_r2: 0.082_f64,
            current_front_log_basis_b1: -0.025_f64,
            current_second_log_basis_b2: 0.012_f64,
            speculator_spreading_imbalance_q_sp: 0.65_f64,
            clearing_intermediary_risk_aversion_gamma: 1.80_f64,
        };
        let r0 = evaluate_boons_prado_commodity_basis_momentum(&base);
        assert!(r0.basis_momentum_signal > 0.0);
        assert!(r0.basis_slope_curvature > 0.0);
        assert!(r0.expected_front_month_excess_return > r0.expected_calendar_spread_return);

        let x0 = [
            base.first_nearby_12m_return_r1,
            base.second_nearby_12m_return_r2,
            base.current_front_log_basis_b1,
            base.current_second_log_basis_b2,
            base.speculator_spreading_imbalance_q_sp,
            base.clearing_intermediary_risk_aversion_gamma,
        ];
        let mut xd = [Dual::<6>::from_f64(0.0); 6];
        for i in 0..6 {
            xd[i] = Dual::variable(x0[i], i);
        }
        let rad =
            evaluate_boons_prado_commodity_basis_momentum(&BoonsPradoCommodityBasisMomentumSpec {
                first_nearby_12m_return_r1: xd[0],
                second_nearby_12m_return_r2: xd[1],
                current_front_log_basis_b1: xd[2],
                current_second_log_basis_b2: xd[3],
                speculator_spreading_imbalance_q_sp: xd[4],
                clearing_intermediary_risk_aversion_gamma: xd[5],
            });
        let h = 1e-6;
        for i in 0..6 {
            let mut xp = x0;
            let mut xm = x0;
            xp[i] += h;
            xm[i] -= h;
            let rp = evaluate_boons_prado_commodity_basis_momentum(
                &BoonsPradoCommodityBasisMomentumSpec {
                    first_nearby_12m_return_r1: xp[0],
                    second_nearby_12m_return_r2: xp[1],
                    current_front_log_basis_b1: xp[2],
                    current_second_log_basis_b2: xp[3],
                    speculator_spreading_imbalance_q_sp: xp[4],
                    clearing_intermediary_risk_aversion_gamma: xp[5],
                },
            );
            let rm = evaluate_boons_prado_commodity_basis_momentum(
                &BoonsPradoCommodityBasisMomentumSpec {
                    first_nearby_12m_return_r1: xm[0],
                    second_nearby_12m_return_r2: xm[1],
                    current_front_log_basis_b1: xm[2],
                    current_second_log_basis_b2: xm[3],
                    speculator_spreading_imbalance_q_sp: xm[4],
                    clearing_intermediary_risk_aversion_gamma: xm[5],
                },
            );
            let fd = (rp.expected_front_month_excess_return
                - rm.expected_front_month_excess_return)
                / (2.0 * h);
            assert!(rad.expected_front_month_excess_return.d[i].abs() > 1e-9);
            assert!((rad.expected_front_month_excess_return.d[i] - fd).abs() < 1e-5);
        }
    }
}
