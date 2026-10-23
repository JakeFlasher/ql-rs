# Changelog

## [0.2.0] — Institutional analytics

- Remediated second-order IFT curvature, AD expiry jets, tail-safe normal tails, exact Heston calibration Hessians, portfolio rate DV01s, AMC envelope AD, commodity VJPs, Dupire jets, and differential-regression scaling.
- Added analytical American option engines, hybrid Heston/Hull–White pricing, rainbow/spread/basket/cliquet APIs, CMS spread pricing, catastrophe bonds, and Black–Litterman/risk-parity allocation.
- Expanded the statistical benchmark harness to 20 suites and the terminal SDK to 18 logical tabs; regenerated the 20-row benchmark reports.

## [0.1.0] — Release

Initial release of the 23-crate `ql-rs` workspace.

### Foundations and contracts

- `ql-time`: dates, day counts, holiday calendars, schedules, and packed business-day operations.
- `ql-ad`: Dual/Dual2 forward AD, reverse arena tape, SIMD scalar support, and IFT primitives.
- `ql-math`: interpolators, distributions, roots, tridiagonal solvers, quadrature, and optimization.
- `ql-curves`, `ql-vol`, `ql-instruments`: multi-curve calibration, volatility surfaces, swaps, bonds, and SoA legs.

### Models and numerical engines

- `ql-models`: Heston/Bates, Hull–White, rough Heston, COS pricing, and calibration.
- `ql-pde`: non-uniform-grid PDE/ADI pricing with multiple splitting schemes.
- `ql-mc`: Philox, Sobol/RQMC, Brownian bridges, SIMD GBM, QE Heston, and two-stage AMC.
- `ql-levy`: Merton, Kou, Variance-Gamma, CGMY, and COS strip pricing.
- `ql-fixed-income` and `ql-rates`: NSS curves, callable/OAS bonds, convertibles, caps/floors, swaptions, G2++, Cheyette, and CMS convexity.

### Asset classes and risk

- `ql-credit`: CDS, hazard curves, Vasicek LHP, ASB CDO recursion, and CDS options.
- `ql-inflation`: seasonal zero-inflation curves and indexed swaps.
- `ql-fx`: Garman–Kohlhagen and sticky-delta Vanna–Volga.
- `ql-exotics`: barriers, Asians, variance swaps, lookbacks, chooser, and quanto options.
- `ql-commodities`: seasonal futures, Schwartz–Smith, Margrabe, Kirk, and spread analytics.
- `ql-scripting`, `ql-xva`, `ql-portfolio`, `ql-dml`: bytecode payoffs, exposure/XVA/SIMM/DIM/MVA, incremental portfolio risk, differential PCA/surrogates, and correlation AAD.

### Delivery and verification

- `ql-benches`: 18 offline benchmark suites with parity assertions and reproducible reports, including RFR SABR and XCCS IFT.
- `ql-terminal`: unified SDK prelude and dependency-free ANSI/Braille multi-asset workstation.

### Extended dependence and products

- `ql-math`: bivariate-normal Genz quadrature and Gaussian/Clayton/Frank/Gumbel/FGM copulas.
- `ql-rates`: backward-looking RFR caplets, Willems SABR smiles, Hull–White RFR analytics, and SOFR futures convexity.
- `ql-fx`: FX forwards, swap points, and multi-currency CSA cross-currency basis IFT curves.
- `ql-exotics`: double barriers, Geske compound options, and complex choosers.
- `ql-inflation`: YoY inflation optionlets and Black/unit-displaced/Bachelier cap-floor models.
- `ql-fixed-income`: amortizing bonds, CPI-linked bonds, and agency MBS prepayment/OAS analytics.
- `ql-portfolio`: coherent Acerbi–Tasche Expected Shortfall and Euler AAD/scenario allocation.

The release includes comprehensive mathematical integration suites and current benchmark output in `knowledge_base/08_BENCHMARK_RESULTS.md`.
