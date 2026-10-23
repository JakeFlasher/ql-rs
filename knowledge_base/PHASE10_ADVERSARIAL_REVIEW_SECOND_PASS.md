# Phase 10 adversarial review — second pass and newly confirmed gaps

This pass re-audits the first mathematical report against the current tree,
searches for placeholders/semantic no-ops, and checks that Phase 10 APIs have
real production workflows rather than only compile-time surface area.

## Newly confirmed critical/high findings

### C1 — Cat-bond Monte Carlo is a deterministic alias, not Monte Carlo

`crates/ql-fixed-income/src/cat_bond.rs:155` implements
`price_mc(&self, _n_paths, _seed)` as `self.price_analytical()` for every peril,
including `PoissonPareto`. This silently returns the exponential analytical
answer for a Pareto model and ignores path count, seed, severity sampling,
discounting and sampling error. It directly violates the Phase 10 contract.

**Required remediation:** implement compound-Poisson event simulation with a
stateless Philox/Sobol stream, sample Poisson counts and Exponential/Pareto
severities, apply Digital/Proportional surviving-notional payoff at each coupon,
return standard errors in the valuation type, and test Pareto MC convergence,
seed reproducibility and lambda→0.

### C2 — Hybrid Heston–Hull–White is an averaging placeholder, not H1-HW COS

`crates/ql-models/src/hybrid.rs:82-151` builds a classical Heston COS price,
builds a BSM price from `theta_v*expiry + hwv`, then returns their arithmetic
average when `eta_r != 0`. `h1_hw_char_func_unit_spot` merely multiplies a Heston
CF by a scalar Gaussian factor and omits the coupled Heston/rate/correlation
terms. This cannot match Grzelak–Oosterlee or QuantLib's analytic H1HW engine,
and it makes `rho_sr` influence only a variance scalar rather than the true
cross terms.

**Required remediation:** implement the coupled affine H1-HW Riccati/characteristic
function and use it in a genuine COS expansion. Add reductions `eta_r→0`,
`xi_v→0`, correlation-sign tests, CF continuity and parameter Dual Greeks versus
central differences/QuantLib.

### C3 — Basket/rainbow Phase 10 formulas are incomplete and several APIs are
semantic no-ops

`ql-exotics/src/basket_rainbow.rs`:

* `TwoAssetCorrelationEngine::price` (lines ~72–96) computes a sum of two
  vanilla-like terms with `(k1+k2)` and returns the same value for calls and
  puts (`let _ = is_call`). It is not Zhang's two-asset correlation option and
  fails put-call parity.
* `DengLiZhouSpreadBasketEngine::price_spread` uses a Kirk-like single
  approximation, not the required Li–Deng–Zhou `J0+J1+J2` boundary expansion;
  its `K=0` branch is an exchange approximation but is not cross-checked against
  Margrabe under all carry conventions.
* `StulzRainbowEngine::price` has an ad-hoc `expected_min` and only one joint
  CDF term; it is not a full Stulz branch implementation for all min/max call
  and put cases.
* `AnalyticCliquetEngine` and `WriterExtensibleOptionEngine` are already
  identified as placeholders in the first report, and this pass confirms their
  `is_call` no-op branches.

**Required remediation:** transcribe the QuantLib branches, preserve option
orientation, use bivariate CDFs with correct transformed correlations and
discounts, and add golden/limit/put-call/Greek tests.

### H1 — American engines are not distinct BAW, BS93, BS2002 and Ju–Zhong
engines

`crates/ql-instruments/src/american.rs` routes all Bjerksund–Stensland and
JuQuadratic calls directly to `BaroneAdesiWhaleyEngine::price`. `critical_spot`
uses a simplified call-only approximation and puts are transformed without a
full put symmetry/parameter check. This means the requested independent engines
and Ju second-order correction are absent despite the API names.

**Required remediation:** implement each referenced approximation separately,
validate free-boundary conditions, add QuantLib/BAW table golden values and
Dual/Dual2 Greek tests. Keep fallback/root failures explicit instead of
`unwrap_or(seed)`.

### H2 — Black–Litterman posterior is mathematically incomplete and f64-only

`crates/ql-portfolio/src/allocation.rs:70-105` returns the prior covariance
unchanged, computes posterior means in f64, and solves an unconstrained weight
vector without the documented fully-invested normalization. The canonical
posterior covariance update and AD/IFT sensitivity lifting requested by Phase 10
are not present. `solve_f64` replaces tiny pivots with `1e-14`, hiding singular or
indefinite covariance matrices.

**Required remediation:** implement the posterior covariance convention
explicitly, validate SPD matrices, expose errors instead of silent pivot repair,
and add canonical He–Litterman two-asset tests plus Dual/Tape sensitivities.

### H3 — Risk-budgeting AD lift does not actually lift covariance/budget jets

`risk_budgeting_weights` converts covariance and budgets to f64 for CCD, then
constructs an approximate correction using frozen primal residuals. It does not
form the requested generic `r_T` from AD covariance/budget inputs and solve the
multivariate IFT Hessian system for all jets. The current `Dual` behavior is
therefore not a valid sensitivity implementation.

**Required remediation:** retain the f64 root only as the primal anchor, form
the full generic residual in `T`, apply one cached inverse Hessian to the generic
residual vector, and test every covariance/budget lane against re-solving CCD.

### M4 — Existing root-lifting and fallback patterns remain production hazards

Second-pass search confirms these independent silent fallbacks remain:

* `ql-rates/src/g2.rs:279` uses `unwrap_or(0.0)` when the nonlinear root solve
  fails.
* `ql-exotics/src/quanto_chooser.rs:87` uses the spot as a failed critical root.
* `ql-instruments/src/american.rs:32` uses the seed when BAW root solving fails.
* `ql-fixed-income/src/callable_bond.rs` and `ql-credit/src/cdo.rs` have solver
  fallbacks/implicit clamps that turn numerical failure into a plausible price.

These should become typed errors or explicit convergence diagnostics; tests must
exercise extreme/ill-conditioned inputs and assert failure rather than accepting
silent values.

### M5 — Phase 10 benchmark integrity still needs an independent oracle

Suites 19/20 now return 20 rows and pass parity, but Suite 20's baseline computes
the requested 90 CCD bumps only as a side effect and returns a scalar sum that is
not a full portfolio allocation checksum; the SOTA returns only `sum(weights)`.
Suite 19 similarly sums heterogeneous price/Greek units and scales both by an
arbitrary `1_000`. These checks can pass while individual risk allocations are
wrong.

**Required remediation:** define a vector checksum (all weights/Greeks), compare
the complete vector with stable units, and independently verify each lane before
timing. Add invariant checks (weights sum to one, RC budgets, American parity).

## Second-pass engineering observations

* `price_mc` placeholder and multiple `let _ = is_call` no-op branches are
  detectable with a simple semantic lint; add a CI grep or deny-list for public
  pricing methods whose arguments are intentionally unused.
* APIs frequently use `unwrap_or` on numerical solvers rather than exposing a
  convergence/error type. This is inconsistent with institutional production
  workflows and should be standardized across crates.
* New Phase 10 tests are mostly smoke/limit tests; the required published-table
  golden values, high-correlation/high-volatility convergence tests, and
  independent QuantLib comparison fixtures are still missing.


## Final delta pass

Cross-checking the math, Rust-engineering and architecture reports produced no new independent issue beyond the findings already listed. The remediation prompt consolidates all unique findings and their regression requirements.
