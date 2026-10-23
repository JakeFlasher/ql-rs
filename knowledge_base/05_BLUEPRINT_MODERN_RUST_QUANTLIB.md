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

`ql-curves` builds OIS discount and IBOR forwarding curves from deposits, FRAs, IRS, and OIS helpers. At a converged pillar vector `theta`, quote sensitivities use Henrard's IFT relation $d\theta/dq=-J_\theta^{-1}J_q$; cached LU factors make repeated VJPs cheap. `ql-vol` wraps Black/Bachelier inversion and SVI/eSSVI/SABR surfaces with Dual/Dual2 sensitivities.

`ql-instruments` resolves swaps, bonds, options, and coupon legs into contiguous arrays (`ResolvedCouponLeg`) for scalar and SIMD valuation. `ql-fixed-income` adds Nelson–Siegel–Svensson zero/forward curves, Black yield-to-price callable bonds, event-aware short-rate trees, OAS duration/convexity, Tsiveriotis–Fernandes convertible-bond induction with cash floors, amortizing sinking-fund bonds, CPI inflation-linked bonds, and PSA/CPR agency MBS pass-through cashflow and prepayment OAS engines. `ql-rates` adds Ibor caplet stripping, physical/cash swaptions, SABR Bachelier pricing, G2++, Cheyette local variance, Linear-TSR CMS convexity, and post-LIBOR backward-looking RFR (SOFR/SONIA/€STR) caplets/floorlets, Willems (2020) effective SABR smiles, Hull–White RFR caplets, and 1M/3M SOFR futures.

### 3.3 Stochastic models and numerical pricing

`ql-models` implements Heston/Bates characteristic functions, COS prices, Hull–White Jamshidian, local volatility, and an Abi Jaber Markovian rough-Heston lift. `ql-levy` implements Merton, Kou, Variance-Gamma, and CGMY characteristic exponents, exact Merton Poisson series, and Fang–Oosterlee COS pricing with Junike–Pankrashkin cumulant truncation; strip pricing reuses characteristic values across strikes.

`ql-pde` solves Heston PDEs with sinh grids, mixed derivatives, boundary conditions, Rannacher damping, and Douglas/Craig–Sneyd/Modified Craig–Sneyd/Hundsdorfer–Verwer ADI sweeps. `ql-mc` provides stateless Philox, Sobol and Brownian bridges, SoA SIMD GBM, Andersen QE Heston, checkpointed pathwise AAD, and two-stage Longstaff–Schwartz AMC. Burley hash-based Owen scrambling and bridge antithetics provide randomized QMC with standard errors.

### 3.4 Asset classes and exotics

`ql-credit` contains ISDA CDS protection/premium/AOD integrals, piecewise hazard bootstrap and IFT Par CS01, Vasicek LHP and Andersen–Sidenius–Basu recursive CDO loss distributions, midpoint tranche spreads, and Black CDS options with front-end protection. `ql-inflation` models multiplicative monthly seasonality, interpolated CPI, ZCIIS/YYIIS, quote IFT deltas, and YoY inflation caps, floors, and collars (Black, Unit-Displaced Black, and Bachelier). `ql-fx` implements Garman–Kohlhagen, sticky-delta strike inversion, Castagna–Mercurio Vanna–Volga Dual3 smiles, FX forwards/swap points, and Constant-Notional Cross-Currency Basis Swaps (`XCCS`) with multi-currency CSA collateral curve bootstrapping and exact IFT cross-currency basis DV01.

`ql-exotics` implements all Reiner–Rubinstein barrier branches with BGK monitoring correction, Ikeda–Kunitomo double-barrier series (`DoubleKnockOut`/`DoubleKnockIn`), Geske compound options (`CallOnCall`, `CallOnPut`, `PutOnCall`, `PutOnPut`) with IFT critical spot, continuous/discrete geometric and Turnbull–Wakeman arithmetic Asians, Demeterfi variance swaps, Goldman–Sosin–Gatto and Conze–Viswanathan lookbacks including the $\lambda \to 0$ L'Hôpital expansion, simple and Rubinstein complex choosers, and quanto options. `ql-commodities` implements seasonal futures IFT deltas, Schwartz–Smith two-factor futures, and Margrabe/Kirk/Spark spreads.

### 3.5 Scripting, XVA, learning, and portfolios

`ql-scripting` lowers declarative products to flat opcodes and evaluates them generically with Fuzzy AAD smoothing for triggers. `ql-xva` simulates `NpvCube`, applies CSA collateral, computes CVA/DVA/FVA, ISDA SIMM, and Chebyshev DIM/MVA profiles. `ql-dml` supplies classical and differential PCA, multi-index differential polynomial regression, and Capriotti–Giles adjoint Cholesky correlation Greeks.

`ql-portfolio` hashes market dependencies, caches calibrated curves and trade results, evaluates dirty trades in Rayon, maps Dual/IFT sensitivities to a portfolio summary, and computes empirical and parametric coherent risk measures (Acerbi–Tasche Expected Shortfall / CVaR, Cornish–Fisher VaR/ES, Markowitz semi-deviation, regret, and Euler homogeneous AAD risk allocation). Together these crates turn individual pricing kernels into scalable scenario and counterparty risk.

### 3.6 Delivery, measurements, and workstation

`ql-benches` is an offline statistical harness with fixed workloads, warmups, medians, IQRs, speedups, and hard parity assertions across 18 suites (`01_calendar_popcnt` through `18_xccs_multicurrency_ift_dv01`). `ql-terminal` re-exports `ql_terminal::prelude::*` and renders a zero-dependency risk workstation: a $2 \times 4$ sub-pixel `BrailleCanvas`, 24-bit ANSI TrueColor ladders, live blotter values, XVA/SIMM/MVA, and calibrated chart series.

## 4. Subsystem Architecture Matrix

| Domain | Crates | Mathematical engines | Verification |
|---|---|---|---|
| Time | `ql-time` | Date serials, calendars, POPCNT, day counts | 01 and time suite |
| AD/math | `ql-ad`, `ql-math` | Dual/Dual2, arena tape, IFT, bivariate normal, copulas, roots/interpolators | AD/math/copula suites |
| Curves/vol | `ql-curves`, `ql-vol` | Multi-curve IFT, Jäckel, SVI/eSSVI/SABR | 03, 04, 10 |
| Rates/fixed income | `ql-instruments`, `ql-fixed-income`, `ql-rates` | SoA legs, NSS, callable/convertible/amortizing/CPI/MBS bonds, G2/Cheyette/CMS/RFR | 02, 17, and rates/fixed-income suites |
| Models/PDE/MC | `ql-models`, `ql-levy`, `ql-pde`, `ql-mc` | COS, ADI, SIMD MC, QE, AMC, RQMC | 05, 06, 14, 15 |
| Asset classes | `ql-credit`, `ql-inflation`, `ql-fx`, `ql-exotics`, `ql-commodities` | CDS/CDO, CPI & YoY caps, FX smiles & XCCS, single/double barriers, compound/chooser/lookback/Asian, futures | 09, 10, 13, 16, 18 and crate suites |
| Risk/automation | `ql-scripting`, `ql-xva`, `ql-dml`, `ql-portfolio` | VM/AAD, XVA/SIMM/DIM, DML, DAG, coherent ES & Euler allocation | 07, 08, 11, 12 and risk suites |
| Delivery | `ql-benches`, `ql-terminal` | Statistical parity, SDK, Braille/ANSI UI | 18-suite harness and terminal suite |
