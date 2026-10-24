# ql-rs Software Architecture and Mathematical Blueprint

## 1. Core tenets

ql-rs uses immutable market snapshots and resolved contracts, generic `Real`/`SimdReal` kernels, IFT at every root/calibration boundary, SoA and arena storage for hot loops, and a content-addressed dependency graph for incremental risk. Pricing functions are explicit `(contract, market, model) -> result` computations; no observer graph or global evaluation date is required.

## 2. Data flow

```text
Date / Calendar / Convention
  -> TradeSpec::resolve -> Resolved SoA cashflows
  -> MarketSnapshot -> calibrated curves and surfaces (IFT Jacobians)
  -> model, PDE, COS, MC, or VM engine (f64 / SIMD / AD)
  -> PortfolioDag / NpvCube / SIMM / terminal output
```

## 3. Crate walkthrough

### 3.1 Foundations: `ql-time`, `ql-ad`, `ql-math`

`ql-time` stores `Date(i32)` serials, civil-date conversions, exchange calendars, day-count conventions, IMM rules, and schedule generation. `BitsetCalendar` packs business-day masks and counts ranges with POPCNT; schedule resolution produces explicit accrual and payment dates.

`ql-ad` defines `Real` and `SimdReal`, so formulas are generic over `f64`, SIMD lanes, `Dual<N>`, `Dual2<N>`, and reverse `TapeVar`. `Dual` carries a fixed gradient, `Dual2` carries a gradient and dense Hessian, and `ArenaTape` stores value/derivative/parent arrays in a rewindable index arena. `ift_lift_1d` and related helpers lift primal solver roots through local Jacobians without recording iterations.

`ql-math` supplies monotone and natural cubic interpolation, univariate and bivariate normal distributions (`bivariate_normal_cdf` via Genz/West/Drezner), elliptical and Archimedean copulas (`GaussianCopula`, `ClaytonCopula`, `FrankCopula`, `GumbelCopula`, `FarlieGumbelMorgensternCopula`), Brent/Newton roots, generic tridiagonal/Thomas solves, Gaussian rules, and Levenberg–Marquardt. These are deliberately scalar-generic and form the numerical substrate for every higher layer.

### 3.2 Curves, volatility, contracts, and rates

`ql-curves` builds OIS discount and IBOR forwarding curves from deposits, FRAs, IRS, and OIS helpers. At a converged pillar vector `theta`, quote sensitivities use Henrard's IFT relation $d\theta/dq=-J_\theta^{-1}J_q$; cached LU factors make repeated VJPs cheap. `ql-vol` wraps Black/Bachelier inversion, SVI/eSSVI/SABR surfaces with Dual/Dual2 sensitivities, Gatheral–Jacquier (2014) Durrleman butterfly-arbitrage density checks (`svi_durrleman_g`), and `CalibratedSabrSmile` / `CalibratedSviSmile` with second-order optimality IFT (`Dual2<N>`) quote sensitivities.

`ql-instruments` resolves swaps, bonds, options, and coupon legs into contiguous arrays (`ResolvedCouponLeg`) for scalar and SIMD valuation, and provides analytical/free-boundary American option approximations (`BaroneAdesiWhaleyEngine`, `BjerksundStensland1993Engine`, `BjerksundStensland2002Engine`, `JuQuadraticEngine`) with `ift_lift_1d` critical boundary and `american_implied_vol` inversion. `ql-fixed-income` adds Nelson–Siegel–Svensson zero/forward curves, Black yield-to-price callable bonds, event-aware short-rate trees, OAS duration/convexity, Tsiveriotis–Fernandes convertible-bond induction with cash floors, amortizing sinking-fund bonds, CPI inflation-linked bonds, PSA/CPR agency MBS pass-through cashflow and prepayment OAS engines, and insurance-linked Catastrophe (`CAT`) bonds (`DigitalNotionalRisk` and `ProportionalNotionalRisk` under Poisson–Exponential, Poisson–Pareto, and compound-Poisson loss processes). `ql-rates` adds Ibor caplet stripping, physical/cash European & Bermudan swaptions (`BermudanSwaption`, `HullWhiteBermudanTreeEngine`, `Fd1dHullWhiteBermudanEngine`), SABR Bachelier pricing, G2++, Cheyette local variance, Linear-TSR CMS convexity, Brigo–Mercurio / Caspers `LognormalCmsSpreadPricer` for CMS spread coupons/caplets/floorlets/collars, and post-LIBOR backward-looking RFR (SOFR/SONIA/€STR) caplets/floorlets, Willems (2020) effective SABR smiles, Hull–White RFR caplets, and 1M/3M SOFR futures.

### 3.3 Stochastic models and numerical pricing

`ql-models` implements Heston/Bates characteristic functions, COS prices, Hull–White Jamshidian, exact `Dual2<2>` Dupire local volatility, an Abi Jaber Markovian rough-Heston lift, and Equity/FX–Rates hybrid engines (`AnalyticBsmHullWhiteEngine` and Grzelak–Oosterlee 2011 `H1-HW` `CosHestonHullWhiteEngine`). `ql-levy` implements Merton, Kou, Variance-Gamma, and CGMY characteristic exponents, exact Merton Poisson series, and Fang–Oosterlee COS pricing with Junike–Pankrashkin cumulant truncation; strip pricing reuses characteristic values across strikes.

`ql-pde` solves Heston PDEs with sinh grids, mixed derivatives, boundary conditions, Rannacher damping, and Douglas/Craig–Sneyd/Modified Craig–Sneyd/Hundsdorfer–Verwer ADI sweeps with full AD time-step sensitivity. `ql-mc` provides stateless Philox, Sobol and Brownian bridges, SoA SIMD GBM, Andersen QE Heston, checkpointed pathwise AAD, and two-stage Longstaff–Schwartz AMC with Glasserman envelope-theorem pathwise AD (`price_with_policy_ad`). Burley hash-based Owen scrambling and bridge antithetics provide randomized QMC with standard errors.

### 3.4 Asset classes and exotics

`ql-credit` contains ISDA CDS protection/premium/AOD integrals, piecewise hazard bootstrap and IFT Par CS01, Vasicek LHP and Andersen–Sidenius–Basu recursive CDO loss distributions, midpoint tranche spreads, and Black CDS options with front-end protection. `ql-inflation` models multiplicative monthly seasonality, interpolated CPI, ZCIIS/YYIIS, quote IFT deltas, and YoY inflation caps, floors, and collars (Black, Unit-Displaced Black, and Bachelier). `ql-fx` implements Garman–Kohlhagen, sticky-delta strike inversion, Castagna–Mercurio Vanna–Volga Dual3 smiles, FX forwards/swap points, and Constant-Notional Cross-Currency Basis Swaps (`XCCS`) with multi-currency CSA collateral curve bootstrapping and exact IFT cross-currency basis DV01.

`ql-exotics` implements all Reiner–Rubinstein barrier branches with BGK monitoring correction, Ikeda–Kunitomo double-barrier series (`DoubleKnockOut`/`DoubleKnockIn`), Geske compound options (`CallOnCall`, `CallOnPut`, `PutOnCall`, `PutOnPut`) with IFT critical spot, continuous/discrete geometric and Turnbull–Wakeman arithmetic Asians, Demeterfi & corridor variance swaps, Broadie–Jain (2008) analytical Heston/Bates variance-of-variance and volatility swap convexity adjustments, Goldman–Sosin–Gatto and Conze–Viswanathan lookbacks including the $\lambda \to 0$ L'Hôpital expansion, simple and Rubinstein complex choosers, quanto options, Stulz (1982) 2-asset Min/Max Rainbow options, Zhang (1995) 2-asset correlation options, Li–Deng–Zhou (2008/2010) 2nd-order boundary & Extended Kirk multi-asset spread/basket options, Levy (1992) moment-matching basket options, Rubinstein forward-start/Cliquet/Performance options, and Longstaff (1990) writer-extensible options. `ql-commodities` implements seasonal futures IFT deltas (single-pass reverse `TapeVar` VJP), Schwartz–Smith two-factor futures, Margrabe/Kirk/Spark spreads, and 1F Ornstein–Uhlenbeck trinomial "Forest of Trees" dynamic programming for commodity swing options (`SwingOptionSpec`) and natural gas storage (`GasStorageSpec`).

### 3.5 Scripting, XVA, learning, and portfolios

`ql-scripting` lowers declarative products to flat opcodes and evaluates them generically with Fuzzy AAD smoothing for triggers. `ql-xva` simulates `NpvCube`, applies CSA collateral, computes CVA/DVA/FVA, Hull–White (2012) / Brigo–Vrins (2018) hazard-rate exposure-coupled Wrong-Way Risk (`compute_xva_with_wwr`), Basel III SA-CCR regulatory Exposure-at-Default (`compute_sa_ccr_ead`), Green–Kenyon–Dennis (2014) Capital Valuation Adjustment (`compute_kva`), ISDA SIMM, and Chebyshev DIM/MVA profiles. `ql-dml` supplies classical and differential PCA, multi-index differential polynomial regression, and Capriotti–Giles adjoint Cholesky correlation Greeks.

`ql-portfolio` hashes market dependencies, caches calibrated curves and trade results, evaluates dirty trades in Rayon, maps Dual/IFT sensitivities to a portfolio summary, computes empirical and parametric coherent risk measures (Acerbi–Tasche Expected Shortfall / CVaR, Cornish–Fisher VaR/ES, Markowitz semi-deviation, regret, and Euler homogeneous AAD risk allocation), and provides institutional portfolio allocation engines (He–Litterman 1999 / Black–Litterman Bayesian posterior, unconstrained and box-constrained Markowitz Mean–Variance & Tangency with active-set KKT implicit differentiation, and Spinu / Bruder–Roncalli Risk Budgeting / Equal Risk Contribution with multivariate IFT Hessian sensitivity lifting). Together these crates turn individual pricing kernels into scalable scenario, counterparty, and asset-allocation risk.

### 3.6 Delivery, measurements, and workstation

`ql-benches` is an offline statistical harness with fixed workloads, warmups, medians, IQRs, speedups, and hard parity assertions across 22 suites (`01_calendar_popcnt` through `22_sabr_svi_ift_and_qp_kkt_aad`). `ql-terminal` re-exports `ql_terminal::prelude::*` and renders a zero-dependency 20-tab risk workstation: a $2 \times 4$ sub-pixel `BrailleCanvas`, 24-bit ANSI TrueColor ladders, live blotter values, XVA/SIMM/MVA/SA-CCR/KVA, and calibrated chart series.

## 4. Subsystem Architecture Matrix

| Domain | Crates | Mathematical engines | Verification |
|---|---|---|---|
| Time | `ql-time` | Date serials, calendars, POPCNT, day counts | 01 and time suite |
| AD/math | `ql-ad`, `ql-math` | Dual/Dual2, arena tape, IFT, bivariate normal, copulas, roots/interpolators | AD/math/copula suites |
| Curves/vol | `ql-curves`, `ql-vol` | Multi-curve IFT & 2nd-order IFT Gamma, Jäckel, SVI/eSSVI/SABR & 2nd-order IFT smile calibration | 03, 04, 10, 22 |
| Rates/fixed income | `ql-instruments`, `ql-fixed-income`, `ql-rates` | SoA legs, American BAW/Bjerksund-Stensland/Ju-Zhong, NSS, callable/convertible/amortizing/CPI/MBS/CAT bonds, Bermudan HW1F Tree/PDE, G2/Cheyette/CMS/CMS-Spread/RFR | 02, 17, 19, 21 and rates/fixed-income suites |
| Models/PDE/MC | `ql-models`, `ql-levy`, `ql-pde`, `ql-mc` | COS, H1-HW & BSM-HW hybrids, ADI, SIMD MC, QE, AMC envelope-AAD, RQMC | 05, 06, 14, 15 |
| Asset classes | `ql-credit`, `ql-inflation`, `ql-fx`, `ql-exotics`, `ql-commodities` | CDS/CDO, CPI & YoY caps, FX smiles & XCCS, single/double barriers, compound/chooser/lookback/Asian/Rainbow/Spread/Cliquet/Extensible, Heston/Bates/corridor variance swaps, futures & Swing/Gas Storage DP | 09, 10, 13, 16, 18 and crate suites |
| Risk/automation | `ql-scripting`, `ql-xva`, `ql-dml`, `ql-portfolio` | VM/AAD, XVA/WWR/SA-CCR/KVA/SIMM/DIM, DML, DAG, coherent ES, Black-Litterman, Box-Markowitz KKT IFT & IFT Risk-Parity | 07, 08, 11, 12, 20, 22 and risk suites |
| Delivery | `ql-benches`, `ql-terminal` | Statistical parity, SDK, 20-tab Braille/ANSI UI | 22-suite harness and terminal suite |
