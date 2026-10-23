# Phase 10 adversarial mathematical review (all workspace crates)

Scope: static review of the 23 crates and Phase-10 additions, with emphasis on quantitative formulas, AD/IFT propagation, numerical conditioning and comparison with standard references (Heston 1993/Little-Heston-Trap, Gatheral 2006, Andersen 2008 QE, Geske 1979, Turnbull--Wakeman, ISDA SIMM and QuantLib implementations). `cargo test --workspace --all-targets` did not expose these issues (most tests use benign parameters and there are no convergence/limit tests for the cases below).

## Findings (ordered by severity)

### CRITICAL

**C1 – Cliquet engine is not pricing a cliquet payoff (wrong payoff and call/put ignored).** `crates/ql-exotics/src/cliquet_extensible.rs:36-46` prices one ATM call (`black_scholes_price(f,f,...)`) and clips that *option price* between `local_floor`/`local_cap`; it never computes `E[min(max(S_{t_i}/S_{t_{i-1}}-1,F),C)]`, and then discards `is_call` (`let _ = is_call`). Floors/caps are on returns, not on the Black price. QuantLib's `CliquetOption` and standard lognormal formulas use a strip of calls/puts at strikes `F(1+F)` and `F(1+C)` plus forwards. Fix with local-return decomposition (or MC), include global floor/cap and sign. Regression: zero vol should equal discounted deterministic clipped returns; `is_call=false` must produce distinct payoff; compare to QuantLib Cliquet.

**C2 – Writer extensible option is a placeholder sum of two unrelated vanillas.** `cliquet_extensible.rs:49-65` returns `BS(S,K1,t1)+BS(S,K2,t2)`. An extensible option has a holder/writer extension decision at `t1`, requiring a max/min continuation value and a bivariate-normal boundary (or tree). Current result violates monotonicity and put-call/extension parity. Replace by Geske-style compound/chooser decomposition and add limits `t2=t1`, `K2=K1`, zero vol.

### HIGH

**H1 – Heston characteristic function has branch discontinuities (no Little-Heston-Trap stabilization).** `crates/ql-models/src/heston.rs:37-48` uses principal complex square root/log with `g=(xi-d)/(xi+d)` directly. As `Re(d)` changes sign, `log(1-g e^{-dT})` jumps, producing oscillatory COS prices/AD Greeks. Gatheral/Lord-Kahl and QuantLib switch to the reciprocal `g` branch (`g=1/G`) and continuous log accumulation. Implement trap form, enforce `Re(d)>=0`, and test CF continuity in `u`, `phi(0)=1`, and calibration Greeks against QuantLib.

**H2 – Heston/Bates COS truncation variance is mathematically incomplete.** `crates/ql-models/src/cos.rs:124-159` sets COS interval width from `mean_var` (integrated expected variance) plus an ad-hoc jump variance. Heston log-price cumulants include vol-of-vol/rho terms in (c_2,c_4); replacing them with (E[\int v,dt]) underestimates tails for high `sigma_v`, long tenor and |rho|≈1. Prices can be negative or converge to the wrong value as `n_terms` grows. Compute analytic cumulants (or cumulant derivatives of the stable CF) and use `a=c1-L√(c2+√c4)`, `b=...`; add convergence tests against semi-analytic Heston integration.

**H3 – Andersen QE Heston branch uses an inconsistent variance/log-price discretization.** `crates/ql-mc/src/paths.rs:230-285` samples `v_{t+dt}` by QE but then uses `v_bar=(v+v')/2` and an extra `rho/sigma_v*(v'-v-kappa(theta-v_bar)dt)` correction. This is neither Andersen's QE-M (which has a specific ψ-dependent log-spot correction) nor exact Broadie--Kaya; it introduces drift bias and can produce wrong forward prices. The `sigma_v=0` branch also uses trapezoidal variance instead of exact integrated CIR (line ~244). Implement QE-M equations from Andersen (2008), clamp/check Feller edge cases, and test martingale `E[S_T]=S_0e^{(r-q)T}` and convergence to Black when `sigma_v=0`.

**H4 – ASB tranche-loss quadrature is renormalized by an ad-hoc first-moment ratio.** `crates/ql-credit/src/cdo.rs:176-187` multiplies the tranche loss by `target_total/quadrature_total`. This changes the nonlinear tranche distribution, can produce loss outside [0,detachment], and is not Andersen-Sidenius-Basu recursion. A quadrature error must be reduced by adaptive/Gauss-Hermite order or exact conditional recursion, not rescaling only the first moment. Remove normalization; add conservation tests (`sum expected issuer losses`), bounds, and comparison with QuantLib GaussianLHPLossModel/ASB.

**H5 – Continuous/discrete cliquet and forward-start discounting are not consistently dated.** In `cliquet_extensible.rs:40-43`, `f=S exp((r-q)t_i)` is used as if a forward at reset but the Black result is additionally discounted from `t_i` via `exp(-r t_{i+1})`; for a local return payoff the conditional value is proportional to (S_{t_i}), requiring discount from `t_i` and a forward ratio, not an ATM call on `f`. Build a date-unit test under flat curves where each reset's expected return is known analytically.

**H6 – SIMM implementation is not ISDA SIMM aggregation.** `crates/ql-xva/src/simm.rs:24-81` computes each bucket `K`, then clips signed sum `S` to `[-K,K]`; ISDA uses prescribed (S_b=\sum WS), (K_b=\sqrt{\sum WS_iWS_j\rho_{ij}}), and inter-bucket (\rho_{bc}\) with sign handling plus risk-class weights, and concentration factor (\sqrt{|s|/\mathrm{threshold}}) only above threshold. No validation of positive-semidefinite correlations or risk-class aggregation is present. Compare to ISDA SIMM v2.6 pseudocode and QuantLib `Simm` tests; add one-/two-bucket golden vectors.

### MEDIUM

**M1 – SIMD normal CDF is only Abramowitz-Stegun erf (~1e-7 absolute error), unsuitable for Greeks.** `crates/ql-ad/src/scalar.rs:228-254,267-296` uses the 5-term `erf` approximation (max error about 1.5e-7), while scalar path uses `libm::erfc`. SIMD pricing therefore disagrees materially in deep tails and finite-difference Greeks. Use an accurate vectorized erfc (e.g., Cody rational approximants, <=1e-15) or document/tolerate the error; add tail tests at |x|=6..10 against libm.

**M2 – AD square-root at zero silently returns derivative 0.** `crates/ql-ad/src/dual.rs:194-207` sets `sqrt'(0)=0` although derivative is singular. This masks invalid local-vol/variance sensitivities and propagates false zero vegas. Return an explicit non-differentiable/error status, or use a documented subgradient only at payoff kinks (not analytic model roots); test `sqrt(eps)` scaling.

**M3 – `Dual::pow`/`Dual2::pow` assumes positive base and divides by base.** `ql-ad/src/dual.rs` (`fn pow`, around lines 230 and 470) computes `ln(self.val)` and `.../self.val` unconditionally. Negative/zero bases (common for shifted rates and SABR extrapolation) yield NaN/inf derivatives even when exponent is an integer. Restrict domain with checked errors or implement integer-power paths; add zero/negative-base tests.

**M4 – IFT root lifting is only a one-step local correction and can amplify solver residuals.** `crates/ql-ad/src/ift.rs:4-10` evaluates residual once and applies Newton plus a scalar `f_xx` correction. If the scalar root tolerance is not near machine precision (implied-vol and compound roots use ~1e-12), AD value and second-order jets retain (O(r/f_x)) error; near vega=0 this explodes. Require residual/conditioning checks, use safeguarded Newton or implicit derivatives evaluated at the converged root, and return a failure status when |f_x| is small. Add deep-ITM/OTM implied-vol AD tests.

**M5 – `inv_norm_cdf` AD tail handling is unstable.** `crates/ql-math/src/distributions.rs:143-171` clamps only in CDO caller; direct `inv_norm_cdf(T)` accepts p≈0/1, computes Newton with `pdf≈0`, and `ift_lift_1d` divides by that pdf. Return ±∞/error or asymptotic tail approximation and test p∈{1e-15,1−1e-15} with derivative bounds.

**M6 – Inflation CPI interpolation denominator is off by one.** `crates/ql-inflation/src/lib.rs:68-70` uses `(day-1)/dim`; interpolation between two lagged CPI observations over `dim` days should reach the second value at day `dim`, i.e. `(day-1)/(dim-1)` (or clearly define endpoint convention). Current value is biased low at month end. Add day=1/day=dim tests and compare QuantLib interpolation.

**M7 – Curve bootstrap helper is not a bootstrap of quoted instruments.** `crates/ql-curves/src/bootstrap.rs:159-174` maps every quote to `D=1/(1+q t)` regardless of quote type, day count, accrual schedule or prior pillars. This is only a simple-compounded zero-rate conversion; calling it `bootstrap_piecewise_curve`/`bootstrap_global_multicurve` can silently miscalibrate deposits/OIS/IRS. Use helper residual solves sequentially (as in later calibrated system) and test par rates reproduce quotes to <1e-12.

**M8 – Callable-bond tree does not enforce risk-neutral probability bounds.** `crates/ql-fixed-income/src/convertible_bond.rs:21-24` computes (p=(e^{(r-q)dt}-d)/(u-d)) without checking `0≤p≤1`; high vol/large dt gives invalid probabilities and arbitrage. Similar recombining tree in `callable_bond.rs` uses fixed 0.5 probabilities regardless of short-rate drift. Use calibrated Ho-Lee/Hull-White transition probabilities (or clamp/refine dt), assert bounds, and add martingale tests.

**M9 – `G2++` swaption integration silently falls back to an incorrect root and clips correlation.** `crates/ql-rates/src/g2.rs:224-264` clips `rho_xy` to ±0.999999 and Newton-refines a quadrature root regardless of monotonicity; fallback `y0=0` when Brent fails can materially change payer/receiver prices. Return a solver error, bracket using monotonic bond value bounds, and test against QuantLib G2 analytic engine for extreme rho/low vol.

**M10 – Double-barrier series has no convergence/error control and clamps negative prices.** `crates/ql-exotics/src/double_barrier.rs:117-149` truncates at user `series_terms` and applies `.max_val(0)`, hiding divergent/oscillatory image sums. Add adaptive tail bound, validate barrier ordering/spot domain, and test in/out parity and convergence as terms increase against Reiner--Rubinstein/QuantLib.

### LOW / robustness

**L1 – `DiscountCurve` permits negative/greater-than-one discounts and extrapolates unconstrained slopes.** `ql-curves/src/curves.rs:51-58,85-103` checks only `d>0`; extreme nodes can imply exploding forwards and negative discount arbitrage. Add monotonicity/upper-bound checks (or explicit negative-rate mode), and extrapolation policy.

**L2 – Bivariate-normal CDF switches approximations at |rho|=0.925 with no continuity/error test.** `ql-math/src/distributions.rs:46-126` uses separate Genz/Drezner branches; derivative/price can jump at threshold. Blend or use one robust algorithm and test continuity in rho and tails.

**L3 – Barrier/discrete BGK adjustment lacks validity checks.** `ql-exotics/src/barrier.rs:37-53` applies BGK shift for any observation count, including one/few observations and barriers at/through spot; BGK is asymptotic for fine monitoring and can worsen prices. Add discrete-monitoring MC benchmarks and disable/flag when `N` small.

## Suggested minimum regression suite

1. Golden values vs QuantLib for Heston (CF/COS), G2 swaption, cliquet, writer-extensible, ASB/LHP CDO and double barriers.
2. Limit tests: zero vol, zero maturity, sigma_v→0 (Black), rho→±1, barriers touching spot, one-point curves.
3. Martingale/parity checks for every MC engine and in/out, call/put and cap/floor parity.
4. AD checks (first and second derivatives) against complex-step/central differences, especially tails and IFT roots with small vegas.
5. Convergence tests in COS/series/quadrature terms and explicit error reporting instead of positivity clamps.

## Reference points

* Heston (1993), “A Closed-Form Solution for Options with Stochastic Volatility”; Lord & Kahl (2006), “Complex Logarithms in Heston-like Models” (Little-Heston-Trap).
* Andersen (2008), “Simple and Efficient Simulation of the Heston Stochastic Volatility Model” (QE/QE-M).
* Gatheral, *The Volatility Surface* (2006), Dupire local-vol denominator.
* Geske (1979), compound options; Rubinstein (1991), chooser/extensible options.
* Andersen–Sidenius–Basu (2003) and ISDA SIMM methodology for structured-credit/SIMM aggregation.

**M11 – Black--Litterman posterior covariance and AD sensitivities are incomplete.** `crates/ql-portfolio/src/allocation.rs:80-128` computes posterior mean in `f64`, returns the *prior* covariance unchanged, and solves weights via a pivot that replaces any negative pivot by `+1e-14` (`solve_f64`). Canonical BL posterior covariance is \(\Sigma+\tau\Sigma P^T(P\tau\Sigma P^T+\Omega)^{-1}P\Sigma\) (or its chosen convention), and all covariance/view inputs should be validated SPD. Preserve AD jets or explicitly document f64-only behavior; add canonical BL two-asset golden test and singular/indefinite matrix errors instead of `unwrap()`.

**M12 – Empirical ES/VaR conventions are easy to misapply and lack input validation.** `ql-portfolio/src/risk_measures.rs:44-75` interprets observations as PnL and returns positive loss, but accepts arbitrary confidence only via assertions and uses the lower-tail index `floor(n*(1-c))` without documenting interpolation at exact integer mass. Add explicit tail convention/types, empty/NaN checks, and Acerbi--Tasche reference vectors.

## Final delta pass

A final cross-report/static delta pass was performed against the Rust and architecture reports and the second-pass semantic search. No additional independent mathematical defect was found beyond the findings recorded above and in `PHASE10_ADVERSARIAL_REVIEW_SECOND_PASS.md`; the remaining items are intentionally consolidated there rather than duplicated.
