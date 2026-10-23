# Literature and References

The 34 papers below are the mathematical and systems references for ql-rs. Each entry names the principal equation or algorithm and the implementing crate/module.

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

The engines also follow Black (1976), Heston (1993), Merton (1976), Kou (2002), Madan–Carr–Chang (1998), Carr–Geman–Madan–Yor (2002), Fang–Oosterlee (2008), Reiner–Rubinstein (1991), Broadie–Glasserman–Kou (1997), Ikeda–Kunitomo (1992), Geske (1979), Rubinstein (1991), Drezner (1978), Genz (2004), Turnbull–Wakeman (1991), Demeterfi–Derman–Kamal–Zou (1999), Goldman–Sosin–Gatto (1979), Conze–Viswanathan (1991), Margrabe (1978), Kirk (1995), Schwartz–Smith (2000), Tsiveriotis–Fernandes (1998), Vasicek (2002), Andersen–Sidenius–Basu (2003), Castagna–Mercurio (2007), Lyashenko–Mercurio (2019), Willems (2020), Rockafellar–Uryasev (2000), Acerbi–Tasche (2002), and ISDA SIMM methodology.
