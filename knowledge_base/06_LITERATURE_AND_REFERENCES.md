# Literature and References

The 59 papers below are the mathematical and systems references for ql-rs. Each entry names the principal equation or algorithm and the implementing crate/module.

## AAD, IFT, DML, and contract languages

1. **`01_Giles_Glasserman_2005_Smoking_Adjoints_Fast_MC_Greeks.pdf`** — AAD for Monte Carlo: reverse propagation makes a path payoff gradient independent of the number of risk inputs. Implemented by ql-ad::ArenaTape, ql-mc pathwise tape, and ql-benches suite 06.

2. **`02_Capriotti_Giles_2010_Fast_Correlation_Greeks_AAD.pdf`** — Adjoint Cholesky factorization and correlation derivatives. Implemented in ql-dml::cholesky_aad and verified by suite 11.

3. **`03_Henrard_2011_AAD_Calibration_Implicit_Function_Theorem.pdf`** — IFT calibration relation dtheta/dq = -Jtheta^-1 Jq. Used by ql-curves, ql-credit, ql-inflation, ql-fx, ql-rates, and ql-commodities.

4. **`04_Huge_Savine_2020_Differential_Machine_Learning.pdf`** — Differential regression augments value observations with pathwise Greeks. Implemented by ql-dml::regression.

5. **`05_Huge_Savine_2025_Axes_That_Matter_Differential_PCA.pdf`** — Differential PCA extracts payoff-relevant gradient axes rather than state covariance axes. Implemented by ql-dml::pca.

6. **`06_Moses_Churavy_2020_Enzyme_LLVM_Autodiff.pdf`** — LLVM-level reverse differentiation motivates scalar-generic kernels and explicit tape boundaries in ql-ad.

7. **`07_Liang_et_al_2025_ad_trait_Rust_Autodiff.pdf`** — Trait-based Rust AD design informs ql-ad::Real, Dual, Dual2, TapeVar, and SimdReal.

8. **`08_Margossian_2019_Review_of_Automatic_Differentiation.pdf`** — Forward/reverse mode complexity and checkpointing provide the AD design criteria used throughout ql-rs.

9. **`14_Bianchetti_2009_Two_Curves_One_Price_MultiCurve.pdf`** — OIS discounting and tenor forwarding separation informs ql-curves and ql-instruments.

10. **`15_PeytonJones_Eber_Seward_2000_Composing_Contracts.pdf`** — Algebraic contract combinators motivate ql-scripting AST construction and event semantics.

11. **`16_Bahr_Berthold_Elsman_2015_Symbolic_Financial_Contracts.pdf`** — Symbolic causality and contract rewriting motivate ql-scripting declarative products.

12. **`19_Huge_Savine_2025_Differential_ML_With_A_Difference.pdf`** — Differential ML with AAD risk labels informs ql-dml and portfolio surrogate workflows.

## PDE, Monte Carlo, COS, and RQMC

1. **`09_Haentjens_intHout_2012_ADI_Heston_Hull_White_PDE.pdf`** — Douglas, Craig–Sneyd, MCS, and HV ADI schemes with mixed derivatives. Implemented by ql-pde::adi.

2. **`10_Duering_et_al_2015_High_Order_ADI_Stochastic_Vol.pdf`** — High-order compact mixed-derivative discretization informs ql-pde non-uniform stencils and cross terms.

3. **`11_Oancea_et_al_2019_Functional_Accelerating_American_MC.pdf`** — Parallel regression and functional-array AMC patterns inform ql-mc::AmcEngine.

4. **`17_Junike_Pankrashkin_2022_Precise_COS_Truncation_Range.pdf`** — Cumulant-based COS truncation range informs ql-models and ql-levy COS engines.

5. **`18_Heston_Simulation_Methods_Comparison_2025.pdf`** — QE versus Euler simulation in Feller-violating regimes informs ql-mc::HestonQeSimulator.

6. **`21_Lakhany_Zhang_2021_ISDA_SIMM_Least_Squares_MC.pdf`** — Least-squares SIMM aggregation informs ql-xva::simm and ql-dml risk aggregation.

7. **`22_Abi_Jaber_2019_Lifting_the_Heston_Model.pdf`** — Exponential-factor Markovian lift of the fractional kernel informs ql-models::rough.

8. **`28_Burley_2020_Practical_Hash_Based_Owen_Scrambling.pdf`** — Laine–Karras hash-based nested Owen scrambling informs ql-mc::BurleyOwenSobol and suite 14.

## XVA, rates, rough volatility, credit, inflation, FX, commodities, RFR, and risk

1. **`12_Kenyon_Green_2014_Efficient_XVA_Management_Regression.pdf`** — Regression exposure cubes unify CVA, DVA, FVA and margin valuation. Implemented in ql-xva.

2. **`13_Albanese_Crepey_2020_XVA_Analysis_Balance_Sheet.pdf`** — Netting-set and balance-sheet XVA architecture informs ql-xva::NpvCube and CSA metrics.

3. **`20_Zeron_Ruiz_2018_Dynamic_Initial_Margin_SIMM_Chebyshev.pdf`** — Chebyshev model-space interpolation for dynamic SIMM/MVA informs ql-xva::chebyshev.

4. **`23_El_Euch_Rosenbaum_2019_Characteristic_Function_Rough_Heston.pdf`** — Fractional Riccati representation and short-time skew inform ql-models rough COS tests.

5. **`24_Kerkhof_2005_Inflation_Derivatives_Explained.pdf`** — CPI lags, daily interpolation, seasonality, ZCIIS, YYIIS, and YoY inflation caps/floors inform ql-inflation.

6. **`25_Bossens_Rayee_Skantzos_Deelstra_2010_Vanna_Volga_FX_Derivatives.pdf`** — Sticky-delta 25D RR/BF inversion and Vanna–Volga correction inform ql-fx.

7. **`26_Brigo_Capponi_2010_Credit_Default_Swaps_Hazard_Rate_Calibration.pdf`** — Piecewise hazard protection and accrual-on-default integrals inform ql-credit.

8. **`27_Gairat_Gorovoy_Shcherbakov_2025_Explicit_Local_Vol_Cheyette_Models.pdf`** — Explicit Cheyette local variance equations inform ql-rates::cheyette.

9. **`29_Fontana_2023_Caplet_Pricing_Affine_Models_Alternative_RFR.pdf`** — Forward-looking vs. backward-looking RFR caplets/floorlets, Jensen's inequality $V_{\text{cpl}}^{\text{bwd}} \ge V_{\text{cpl}}^{\text{fwd}}$, term-basis caplets, and 1M/3M RFR futures inform ql-rates::rfr.

10. **`30_Acerbi_Tasche_2002_On_the_Coherence_of_Expected_Shortfall.pdf`** — Universal coherence of Expected Shortfall on general (continuous and discontinuous) loss distributions and exact quantile tail correction inform ql-portfolio::risk_measures.

11. **`31_Baran_2016_Post_Crisis_Valuation_OIS_Cross_Currency_Basis_Swaps.pdf`** — Collateralized OIS discounting, Constant-Notional Cross-Currency Basis Swap (`XCCS`) curve bootstrapping ($DF_{f,\text{CSA}}$ from $DF_{d,\text{OIS}}$, $DF_{f,\text{OIS}}$, and $BS_n$), and FX swap points inform ql-fx::xccs and suite 18.

12. **`32_Acerbi_Tasche_2002_Expected_Shortfall_Natural_Coherent_Alternative.pdf`** — Order-statistic estimator, sub-additivity, and Euler component risk allocation for Expected Shortfall (CVaR) inform ql-portfolio::risk_measures.

13. **`33_Willems_2020_SABR_Smiles_for_RFR_Caplets.pdf`** — Closed-form effective SABR parameters $(\hat{\alpha}, \beta, \hat{\rho}, \hat{\nu})$ for backward-looking compounded-in-arrears RFR caplets both before and inside the accrual period inform ql-rates::rfr and suite 17.

14. **`34_Turfus_RomeroBermudez_2023_Analytic_RFR_Option_Pricing_Smile_Skew.pdf`** — Hull–White and quadratic-skew short-rate pricing kernels and effective term variance for compounded RFR caplets and swaptions inform ql-rates::rfr.

15. **`35_Grzelak_Oosterlee_2011_Heston_Stochastic_Interest_Rates.pdf`** — Closed-form affine approximation of the Heston–Hull–White (`H1-HW` and `H2-HW`) characteristic function with non-zero equity–rate correlation $\rho_{Sr}$ informs `ql-models::hybrid::CosHestonHullWhiteEngine`.

16. **`36_Li_Deng_Zhou_2008_Multi_Asset_Spread_Option_Pricing.pdf`** — Second-order exercise-boundary approximation and Extended Kirk formula for 2-asset and multi-asset spread and basket options inform `ql-exotics::basket_rainbow::DengLiZhouSpreadBasketEngine`.

17. **`37_Choi_2018_Sum_of_All_BSM_Models_Basket_Spread_Options.pdf`** — Factor-rotated Gauss–Hermite quadrature and moment-matching benchmarks for multi-asset BSM basket and spread options inform `ql-exotics::basket_rainbow`.

18. **`38_GriveauBillion_Richard_Roncalli_2013_Fast_Risk_Parity_Portfolios.pdf`** — Cyclical Coordinate Descent (CCD) and Spinu Newton–Cholesky optimization on the strictly convex risk-budgeting potential $\frac{1}{2} y^\top \Sigma y - \sum_i b_i \ln y_i$ inform `ql-portfolio::allocation` and suite 20.

19. **`39_Bruder_Roncalli_2012_Risk_Budgeting_Approach.pdf`** — Risk Budgeting Euler decomposition $RC_i = w_i (\Sigma w)_i / \sqrt{w^\top \Sigma w}$, closed-form bivariate solution, and volatility ordering $\sigma_{\text{MV}} \le \sigma_{\text{RB}} \le \sigma_{\text{WB}}$ inform `ql-portfolio::allocation`.

20. **`40_He_Litterman_1999_Intuition_Behind_Black_Litterman_Portfolios.pdf`** — CAPM equilibrium implied excess returns $\Pi = \delta \Sigma w_{\text{mkt}}$ and Bayesian posterior view-tilting $(\mu_{\text{BL}}, \Sigma_{\text{BL}})$ inform `ql-portfolio::allocation::BlackLittermanSpec`.

21. **`41_BaroneAdesi_Whaley_1987_Efficient_Analytic_Approximation_American_Options.pdf`** — Quadratic free-boundary early-exercise approximation and critical spot $S^*$ root-finding for American calls and puts inform `ql-instruments::american` and suite 19.

22. **`42_Miemiec_2014_Pricing_CMS_Spread_Options_QuantLib.pdf`** — Shifted-lognormal Gauss–Hermite conditioning and Linear-TSR convexity adjustment for CMS spread coupons, caplets, floorlets, and collars inform `ql-rates::cms_spread`.

23. **`43_Blondel_et_al_2021_Efficient_Modular_Implicit_Differentiation.pdf`** — Second-order optimality Implicit Function Theorem $\partial \theta^* / \partial q = -(\nabla^2_{\theta\theta} L)^{-1} \nabla^2_{\theta q} L$ including the residual-Hessian term $\sum_i r_i \nabla^2_{\theta\theta} \hat{y}_i(\theta^*)$ and active-set KKT differentiation inform `ql-models::calibration`, `ql-vol::surface` (`CalibratedSabrSmile`, `CalibratedSviSmile`), and `ql-portfolio::allocation` (`BoxConstrainedMarkowitzSpec`).

24. **`44_Jacquier_et_al_2024_Deep_Calibration_Rough_Vol.pdf`** — Rough volatility calibration regularizers and Markovian kernel parameterizations inform `ql-models::rough`.

25. **`45_Agrawal_et_al_2019_Differentiable_Convex_Optimization_Layers.pdf`** — Differentiating through constrained quadratic/convex optimization layers via active-set KKT conditions informs `ql-portfolio::allocation::try_box_constrained_markowitz_weights` and suite 22.

26. **`46_Gatheral_Jacquier_2014_Arbitrage_Free_SVI_Volatility_Surfaces.pdf`** — Durrleman butterfly arbitrage condition $g(k) = (1 - \frac{k w'(k)}{2w(k)})^2 - \frac{w'(k)^2}{4}(\frac{1}{w(k)} + \frac{1}{4}) + \frac{w''(k)}{2} \ge 0$, Roger Lee wing slope bound $b(1 + |\rho|) \le 2$, and arbitrage-free SVI calibration inform `ql-vol::surface` (`svi_durrleman_g`, `CalibratedSviSmile`).

27. **`47_Green_Kenyon_Dennis_2014_KVA_Capital_Valuation_Adjustment.pdf`** — Counterparty credit risk capital profiles $K_{\text{CCR}}(t)$, Basel III SA-CCR Exposure-at-Default $\text{EAD} = 1.4(\text{RC} + \text{multiplier}\times\text{AddOn}^{\text{agg}})$, and Capital Valuation Adjustment (`compute_sa_ccr_ead`, `compute_kva`) inform `ql-xva`.

28. **`48_Brigo_Vrins_2018_Disentangling_Wrong_Way_Risk_CVA.pdf`** — Wrong-Way Risk CVA/DVA exposure-intensity coupling and survival-calibrated pathwise hazard rates $\lambda_c(t_k, \omega) = \exp(a_k + b_{\text{wwr}} w(t_k, \omega))$ inform `ql-xva::compute_xva_with_wwr`.

29. **`49_Warin_2012_Hedging_Swing_Contract_on_Gas_Market.pdf`** — 1F mean-reverting Ornstein–Uhlenbeck log-spot trinomial tree calibrated to `CommodityFuturesCurve` and 2D backward dynamic programming ("Forest of Trees") with mandatory take-or-pay exercise constraints inform `ql-commodities` (`SwingOptionSpec`, `price_swing_option`).

30. **`50_Yamakami_Takeuchi_2022_Bermudan_Swaption_Two_Factor_Hull_White.pdf`** — Single-curve and dual-curve (OIS discount / IBOR-RFR forward) Bermudan swaption backward induction and affine bond reconstruction inform `ql-rates::swaption` (`BermudanSwaption`, `HullWhiteBermudanTreeEngine`, `Fd1dHullWhiteBermudanEngine`) and suite 21.

31. **`51_Baeuerle_Riess_2015_Gas_Storage_Valuation_Multinomial_Tree_LSM.pdf`** — Markov Decision Process backward induction over inventory grid $v \in [0, V_{\max}]$ with injection/withdrawal rate limits, proportional loss factors, and unit costs on a recombining spot tree informs `ql-commodities` (`GasStorageSpec`, `price_gas_storage`).

32. **`52_McCrickerd_Pakkanen_2018_Turbocharging_MC_Rough_Bergomi.pdf`** — Hybrid scheme for Riemann–Liouville fractional Brownian motion $\widetilde{W}_t^H = \sqrt{2H}\int_0^t (t-s)^{H-1/2}dW_s^{(1)}$ with exact singular $2\times 2$ head covariance, optimal power-law weights $b_m$, and Romano–Touzi conditional Black–Scholes variance reduction informs `ql-models::rough` (`RoughBergomiParams`, `rbergomi_hybrid_mc_price`) and suite 23.

33. **`53_Cozma_Mariapragassam_Reisinger_2018_SLV_Particle_Method.pdf`** — McKean–Vlasov kernel-smoothing particle calibration of the Gyöngy leverage surface $L(t, S) = \sigma_{\text{Dup}}(t, S)/\sqrt{\mathbb{E}[V_t \mid S_t = S]}$ in hybrid local-stochastic volatility models informs `ql-models::slv` (`calibrate_slv_particle_leverage`, `price_slv_european_mc`).

34. **`54_Reisinger_Tsianni_2025_Particle_System_Heston_SLV.pdf`** — Propagation of chaos, Nadaraya–Watson quartic biweight / Gaussian kernel regularization, and martingale projection for calibrated Heston-type SLV particle systems inform `ql-models::slv`.

35. **`55_Alfonsi_Fruth_Schied_2010_Optimal_Execution_Transient_Impact.pdf`** — Obizhaeva–Wang transient limit-order-book resilience $G(t-s) = \theta_0 e^{-\rho(t-s)}$, U-shaped boundary-block optimal execution schedules, and Almgren–Chriss permanent/temporary impact trajectories inform `ql-portfolio::execution` (`try_almgren_chriss_execution`, `try_obizhaeva_wang_execution`) and suite 24.

36. **`56_Bichuch_Capponi_Sturm_2018_Arbitrage_Free_XVA_Funding.pdf`** — Arbitrage-free nonlinear XVA BSDEs/PDEs under asymmetric borrowing/lending funding rates ($r_f^+ \ne r_f^-$) and collateral rates ($r_c$), decomposing bilateral valuation into CVA, DVA, FCA, and FBA via fixed-point iteration, inform `ql-xva` (`BilateralFundingSpec`, `compute_bilateral_xva_fixed_point`).

37. **`57_Rockafellar_Uryasev_2000_Optimization_of_CVaR.pdf`** — Auxiliary convex functional $F_\alpha(w, \zeta) = \zeta + \frac{1}{(1-\alpha)S}\sum_{s=1}^S [-R_s^\top w - \zeta]_+$ for simultaneous Value-at-Risk ($\zeta^*$) and Expected Shortfall ($\text{CVaR}_\alpha$) minimization informs `ql-portfolio::allocation` (`MeanCvarSpec`, `try_mean_cvar_weights`) and suite 24.

38. **`58_LopezSalas_2024_PDEs_Generalized_Forward_Market_Model.pdf`** — Generalized Forward Market Model (FMM / BGM) for forward-looking IBOR and backward-looking compounded-in-arrears RFR rates with decaying in-accrual volatility $\gamma_i(t) = (T_{i+1}-t)/(T_{i+1}-T_i)$, Hunter–Jäckel–Joshi predictor–corrector drift, and TARN/Callable Range Accrual pricing inform `ql-rates::lmm`.

39. **`59_Tasche_2018_Capital_Allocation_FRTB_Expected_Shortfall.pdf`** — Basel III / IV Fundamental Review of the Trading Book (FRTB MAR20–MAR33) liquidity-horizon-scaled Expected Shortfall $ES = \sqrt{ES_1^2 + \sum_{j \ge 2}(ES_j\sqrt{(LH_j - LH_{j-1})/T_{\text{base}}})^2}$, Sensitivities-Based Method (SBM Delta/Vega/Curvature under Low/Medium/High correlation scenarios + DRC + RRAO), and Euler AAD capital allocation inform `ql-xva::frtb` and suite 24.

## Reference repositories

| Repository | Role in ql-rs design |
|---|---|
| `repos/QuantLib` | Canonical C++ contracts, conventions, analytic engines, and regression references. |
| `repos/RustQuant` | Comparative Rust audit and identification of completeness gaps. |
| `repos/quantsupport` | Unified scalar traits, AAD tape, scripting, and immutable stores. |
| `repos/rateslib` | Manifold-shared Dual/Dual2 and global multi-curve IFT patterns. |
| `repos/stochastic-rs` | SIMD/GPU stochastic-process and PDE patterns. |
| `repos/OptionStratLib` | Rust option analytics and risk API comparisons. |
| `repos/volsurf` | Volatility surface conventions and arbitrage checks. |
| `repos/implied-vol` | Robust Jäckel-style implied-volatility inversion. |
| `repos/CompFinance` | Savine block-arena AAD and path checkpointing. |
| `repos/Scripting` | Flat opcode compilation and Fuzzy AAD payoff evaluation. |

## Classic analytical references

The engines also follow Black (1976), Heston (1993), Merton (1976), Kou (2002), Madan–Carr–Chang (1998), Carr–Geman–Madan–Yor (2002), Fang–Oosterlee (2008), Barone-Adesi–Whaley (1987), Bjerksund–Stensland (1993, 2002), Ju–Zhong (1999), Reiner–Rubinstein (1991), Broadie–Glasserman–Kou (1997), Ikeda–Kunitomo (1992), Geske (1979), Rubinstein (1991), Stulz (1982), Zhang (1995), Levy (1992), Longstaff (1990), Drezner (1978), Genz (2004), Turnbull–Wakeman (1991), Demeterfi–Derman–Kamal–Zou (1999), Goldman–Sosin–Gatto (1979), Conze–Viswanathan (1991), Margrabe (1978), Kirk (1995), Li–Deng–Zhou (2008), Schwartz–Smith (2000), Tsiveriotis–Fernandes (1998), Vasicek (2002), Andersen–Sidenius–Basu (2003), Castagna–Mercurio (2007), Brigo–Mercurio (2006), Lyashenko–Mercurio (2019), Willems (2020), He–Litterman (1999), Bruder–Roncalli (2012), Spinu (2013), Rockafellar–Uryasev (2000), Acerbi–Tasche (2002), and ISDA SIMM methodology.
