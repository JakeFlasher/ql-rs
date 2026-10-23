# 01. Anatomy of C++ QuantLib (`lballabio/QuantLib`): Brilliance & 25-Year Architectural Debt

> **Repository Inspected**: [`/home/jakeshea/rust_quantlib/repos/QuantLib`](file:///home/jakeshea/rust_quantlib/repos/QuantLib) (`HEAD` commit `bf35b19a`, Oct 2026)  
> **Codebase Scale**: **443,754 LOC** across 2,438 `.hpp`/`.cpp` files in `ql/`, plus **134,959 LOC** in `test-suite/` (**578,713 total C++ LOC**).

---

## 1. Quantitative Scale & Module Breakdown of C++ QuantLib

Before asking whether a modern Rust successor is needed, one must appreciate the sheer domain surface area accumulated in C++ QuantLib over 26 years (since 2000). A line-of-code audit of `/home/jakeshea/rust_quantlib/repos/QuantLib/ql` reveals where the complexity actually lives:

| Subsystem (`ql/`) | C++ LOC | Key Responsibilities & Institutional Domain Knowledge |
| :--- | ---: | :--- |
| `ql/math/` | **151,144** | 1D/2D interpolations (cubic splines, monotonic, Akima, SABR, convex monotone), solvers (Brent, Newton, Ridder), optimizers (Levenberg-Marquardt, BFGS, Simplex, Differential Evolution), matrix decompositions (SVD, symmetric Schur, Cholesky, QR), copulas, Sobol/Halton low-discrepancy sequences (Joe-Kuo direction numbers), Gaussian quadratures. |
| `ql/experimental/` | **68,331** | Credit/CDS extensions, CVA/XVA prototypes, exotic barrier/variance/commodity/catastrophe instruments, SABR/ZABR extensions, Gaussian 1D/Markov functional models. |
| `ql/pricingengines/` | **47,118** | Analytic, FD (PDE), Monte Carlo, lattice, and Fourier/Cosine/FFT engines across vanilla, barrier, Asian, basket, cliquet, swaption, cap/floor, bond, credit, and inflation instruments. |
| `ql/termstructures/` | **33,940** | `YieldTermStructure`, `PiecewiseYieldCurve`, `IterativeBootstrap`, `GlobalBootstrap`, `MultiCurve` (cyclic multi-curve solver), Black/Local/Heston/SABR volatility surfaces, swaption vol cubes, caplet stripping, default probability & inflation curves. |
| `ql/models/` | **33,564** | Short-rate (Hull-White, Black-Karasinski, G2++, Cox-Ingersoll-Ross, Generalized Hull-White), Libor Market Model (LMM), Heston, Bates,GJRGARCH, Markov Functional, Gaussian1D/GSR, calibration helpers. |
| `ql/instruments/` | **23,588** | Fixed-income (`Bond`, `FixedVsFloatingSwap`, `VanillaSwap`, `OvernightIndexedSwap`, `CrossCurrencyBasisSwap`, `Swaption`, `CapFloor`, `ForwardRateAgreement`), equity/FX options, credit default swaps, inflation swaps. |
| `ql/methods/` | **20,901** | `finitedifferences/` (1D/2D/3D PDE operators, mesher, ADI schemes: Douglas, Craig-Sneyd, Hundsdorfer-Verwer, TrBDF2, Method of Lines), `montecarlo/` (path generators, Longstaff-Schwartz AMC, Brownian bridge), `lattices/` (binomial/trinomial trees). |
| `ql/time/` | **20,743** | `Date` (Excel-style serial numbers), `Period`, `Frequency`, `Calendar` (45+ national/exchange holiday rules), `DayCounter` (Act/360, Act/365F, Act/Act ISDA/ICMA/AFB, 30/360 US/Eurobond/ISDA, Business252), `Schedule` (stub handling, EOM, CDS/IMM rules). |
| `ql/cashflows/` | **15,921** | `Coupon`, `FixedRateCoupon`, `IborCoupon`, `OvernightIndexedCoupon` (daily compounding/averaging with lockout, lookback, observation shift), `CmsCoupon`, coupon pricers (conundrum/convexity adjustments, Bachelier/Black), `CashFlows` (NPV, BPS, IRR, Z-spread, duration, convexity). |
| `ql/indexes/` | **7,213** | `InterestRateIndex`, `IborIndex`, `OvernightIndex` (SOFR, ESTR, SONIA, TONAR), `SwapIndex`, `InflationIndex`, historical fixing repositories (`IndexManager`). |
| `ql/processes/` | **5,387** | `StochasticProcess`, `GeneralizedBlackScholesProcess`, `HestonProcess`, `HullWhiteProcess`, `GsrProcess`, `Merton76Process`, `BatesProcess`, `StochasticProcessArray`. |
| `ql/currencies/` | **2,840** | ISO 4217 currencies, `ExchangeRate`, `ExchangeRateManager` (triangular FX lookup). |
| `ql/patterns/` + Core | **~13,000** | `Observable`/`Observer`, `LazyObject`, `Singleton`, `Visitor`, `Handle`/`RelinkableHandle`, `Instrument`, `PricingEngine`, `Quote`, `Settings`. |

---

## 2. What C++ QuantLib Got Brilliantly Right (Must Preserve)

Any effort to build a Rust-based successor that ignores **why** C++ QuantLib became the industry reference will fail. QuantLib's enduring strengths are:

1. **Decoupling Contract Definition (`Instrument`) from Numerical Valuation (`PricingEngine`)**:
   - In [`ql/instrument.hpp`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/instrument.hpp) and [`ql/pricingengine.hpp`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/pricingengine.hpp), a `VanillaOption` or `VanillaSwap` contains only contractual terms (schedules, strikes, exercise rules, day-count conventions).
   - It does **not** hardcode whether it is priced by Black-Scholes, Heston FFT, a 2D ADI finite-difference grid, or Sobol Monte Carlo, nor does it hold market curves directly. This separation is fundamental to institutional software engineering.
2. **Exhaustive Market Conventions & Edge-Case Fidelity (`ql/time/`, `ql/cashflows/`)**:
   - Real-world fixed income is governed by intricate conventions: `Actual/Actual (ICMA)` requiring reference period start/end dates ([`ql/time/daycounter.hpp:87-89`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/time/daycounter.hpp#L87-L89)), `OvernightIndexedCoupon` daily compounding with lookback/lockout/observation shifts, end-of-month schedule generation with short/long front/back stubs, and ex-coupon bond dates ([`ql/cashflows/coupon.hpp:43-50`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/cashflows/coupon.hpp#L43-L50)).
3. **Composable Multi-Dimensional Finite Difference Framework (`ql/methods/finitedifferences/`)**:
   - Klaus Spanderen's redesign of QuantLib's FDM framework ([`ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/methods/finitedifferences/solvers/fdmbackwardsolver.hpp)) cleanly separates spatial meshers (`FdmMesher`), differential operators (`FdmLinearOpComposite`), boundary/step conditions (`FdmStepCondition`), and ADI splitting schemes (`Douglas`, `CraigSneyd`, `ModifiedCraigSneyd`, `Hundsdorfer`, `TrBDF2`).
4. **135,000 Lines of Battle-Tested Regression Tests (`test-suite/`)**:
   - Over two decades of institutional bug reports, edge-case calendar rolls, negative interest rate fixes (Bachelier/shifted lognormal), and OIS multi-curve transition corner cases are encoded in `test-suite/`.

---

## 3. The Six Structural Sins of C++ QuantLib (Why It Cannot Be Modernized In-Place)

Designed in 2000 (pre-C++11, pre-multi-core ubiquity, pre-SIMD, pre-AAD, pre-2008 multi-curve crisis), C++ QuantLib is built around 1994 *Gang of Four* stateful Object-Oriented Programming patterns. Inspecting the headers in `repos/QuantLib/ql/` exposes six foundational architectural flaws that cannot be fixed without breaking every class in the library:

### Flaw 1: The Stateful `Observer` / `Observable` / `LazyObject` Pointer Graph

Look at [`ql/patterns/observable.hpp`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/patterns/observable.hpp) and [`ql/patterns/lazyobject.hpp`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/patterns/lazyobject.hpp):

```cpp
// ql/patterns/observable.hpp:63-85
class Observable {
  private:
    typedef std::set<Observer*> set_type;
    set_type observers_;
};

// ql/patterns/observable.hpp:119-159
class Observer {
  private:
    typedef std::set<ext::shared_ptr<Observable>> set_type;
    set_type observables_;
};
```

- **Ubiquitous `virtual` Inheritance & Red-Black Tree Overhead**: Every `Quote`, `TermStructure`, `Index`, `CashFlow`, `Coupon`, `RateHelper`, `PricingEngine`, and `Instrument` inherits via `public virtual Observable, public virtual Observer`. Every single instance owns two `std::set` red-black trees on the heap (`observers_` and `observables_`).
- **Catastrophic Thread-Safety Tax**: In default mode (`!QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN`), QuantLib is **completely thread-unsafe**—even reading `NPV()` on two instruments sharing a yield curve or index on two threads causes data races on `std::set` and `mutable bool calculated_`. When `QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN` is enabled ([`ql/patterns/observable.hpp:290-406`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/patterns/observable.hpp#L290-L406)), every `Observer` allocates a `shared_ptr<Proxy>` with a `std::recursive_mutex` plus its own `std::recursive_mutex`, and every `Observable` adds another `std::recursive_mutex`!
- **Notification Storms**: Updating `Settings::instance().evaluationDate()` or bumping 50 market quotes in a curve triggers cascading virtual `update()` calls through tens of thousands of nodes, acquiring and releasing recursive mutexes at every hop.
- **Cyclic Dependency Fragility**: When Eugene Toder and Peter Caspers added `MultiCurve` in 2025 ([`ql/termstructures/multicurve.hpp:41-84`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/termstructures/multicurve.hpp#L41-L84)) to support cyclic OIS/FX/IBOR multi-curve bootstrapping, they had to introduce complex manual ceremony—separating "internal" `RelinkableHandle`s (using `null_deleter` and disabling `registerAsObserver` to prevent `shared_ptr` memory leaks and infinite `Observer::update()` loops) from "external" `Handle`s.

### Flaw 2: Fake Immutability (`const` Methods Mutating `mutable` Fields)

Look at [`ql/instrument.hpp:107-110`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/instrument.hpp#L107-L110), [`ql/patterns/lazyobject.hpp:131`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/patterns/lazyobject.hpp#L131), [`ql/termstructure.hpp:105-108`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/termstructure.hpp#L105-L108), and [`ql/pricingengine.hpp:72-73`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/pricingengine.hpp#L72-L73):

```cpp
// ql/instrument.hpp:107-109
mutable Real NPV_, errorEstimate_;
mutable Date valuationDate_;
mutable std::map<std::string, std::any> additionalResults_;

// ql/pricingengine.hpp:72-73
mutable ArgumentsType arguments_;
mutable ResultsType results_;
```

- Calling `instrument.NPV() const` or `curve.discount(t) const` appears read-only in C++ signatures, but under the hood it mutates `calculated_`, triggers `bootstrap_.calculate()`, mutates `engine_->arguments_` and `engine_->results_`, and overwrites `NPV_`.
- Consequently, **not even `PricingEngine` or `YieldTermStructure` is stateless or re-entrant**. Two threads cannot price two different swaps using the same `DiscountingSwapEngine` instance, because both threads will overwrite `engine_->arguments_` and `engine_->results_` simultaneously ([`ql/instrument.hpp:147-154`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/instrument.hpp#L147-L154)).

### Flaw 3: Global Mutable Singleton State (`Settings::instance()`)

Look at [`ql/settings.hpp:37-118`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/settings.hpp#L37-L118) and `IndexManager::instance()`:
- The evaluation date (`Settings::instance().evaluationDate()`) and historical index fixings (`IndexManager::instance()`) are global singletons (or per-thread singletons with compile-time session flags).
- In modern risk workflows—such as computing Theta, Horizon PnL, or simulating 100 future exposure dates across 10,000 Monte Carlo paths for XVA—the evaluation date is a **first-class input parameter of the valuation context**, not a process-wide global variable!

### Flaw 4: Global `typedef double Real` Blocks Native AAD & SIMD Coexistence

Look at [`ql/types.hpp:50-82`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/types.hpp#L50-L82):

```cpp
typedef QL_REAL Real;
typedef Real Decimal;
typedef Real Time;
typedef Real DiscountFactor;
typedef Real Rate;
typedef Real Spread;
typedef Real Volatility;
typedef Real Probability;
```

1. **Zero Dimensional Type Safety**: Because `typedef` creates aliases rather than distinct types, passing a `Volatility` into a function expecting a `Rate` or `DiscountFactor`, or passing `Time` where `Rate` is expected, compiles silently without warning.
2. **Global Monomorphic `Real`**: Because classes like `Instrument`, `TermStructure`, `CashFlow`, and `PricingEngine` are **virtual class hierarchies** (and C++ forbids `virtual` member templates!), the scalar type `Real` must be fixed once at compile time for the entire library (`#define QL_REAL double`).
3. **The AAD Dilemma**: In production quantitative finance, you do **not** want every calculation to run on a slow operator-overloading reverse-mode tape! You want:
   - Standard `f64` for untaped root-finding/bootstrapping and primal checks,
   - SIMD vectors (`f64x4` / `f64x8`) for vectorized Monte Carlo and PDE sweeps,
   - Stack-allocated forward dual numbers (`Dual<f64, N>`) for low-dimensional analytic Greeks and curve Jacobians $\frac{\partial H}{\partial z}$,
   - Reverse-mode AAD (`TapeVar` or LLVM Enzyme `#[autodiff_reverse]`) for high-dimensional pathwise sensitivities.
   In C++ QuantLib, because `Real` is a single global typedef behind `virtual` interfaces, you cannot instantiate `BlackCalculator<Dual<double>>` or `HestonProcess<f64x8>` alongside `f64` in the same binary.

### Flaw 5: Runtime Type Erasure (`dynamic_cast` & `std::any`) Instead of Algebraic Types

Look at [`ql/instrument.hpp:147-166`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/instrument.hpp#L147-L166):

```cpp
inline void Instrument::performCalculations() const {
    QL_REQUIRE(engine_, "null pricing engine");
    engine_->reset();
    setupArguments(engine_->getArguments());
    engine_->getArguments()->validate();
    engine_->calculate();
    fetchResults(engine_->getResults());
}

inline void Instrument::fetchResults(const PricingEngine::results* r) const {
    const auto* results = dynamic_cast<const Instrument::results*>(r);
    QL_ENSURE(results != nullptr, "no results returned from pricing engine");
    ...
    additionalResults_ = results->additionalResults; // copies std::map<std::string, std::any>!
}
```

- Nothing at compile time prevents attaching a `BlackSwaptionEngine` to a `VanillaOption`: `instrument.setPricingEngine(engine)` accepts any `ext::shared_ptr<PricingEngine>` and only fails at runtime when `dynamic_cast<VanillaOption::arguments*>(engine_->getArguments())` returns `nullptr`!
- Every valuation copies `std::map<std::string, std::any> additionalResults_`.
- Even `DayCounter` and `Calendar` ([`ql/time/daycounter.hpp:61`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/time/daycounter.hpp#L61), [`ql/time/calendar.hpp:72`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/time/calendar.hpp#L72)) wrap a heap-allocated `ext::shared_ptr<Impl>` and compare equality by allocating two `std::string`s and comparing `d1.name() == d2.name()` ([`ql/time/daycounter.hpp:137`](file:///home/jakeshea/rust_quantlib/repos/QuantLib/ql/time/daycounter.hpp#L137))!

### Flaw 6: Array-of-Structures (AoS) Heap Fragmentation & Pointer Chasing

- A `Leg` in QuantLib is `typedef std::vector<ext::shared_ptr<CashFlow>> Leg;`. Pricing a portfolio of 50,000 30-year swaps means traversing millions of individually heap-allocated polymorphic `CashFlow` / `IborCoupon` objects via `shared_ptr` indirection and vtable dispatch (`cf->amount()`, `cf->date()`, `cf->hasOccurred()`).
- In a Monte Carlo or PDE loop, every call to `Payoff::operator()(Real s)` or `StochasticProcess::evolve(...)` goes through a vtable pointer, preventing the compiler from inlining the payoff/process into the time-stepping loop and completely blocking AVX-512/NEON auto-vectorization or GPU offloading.
