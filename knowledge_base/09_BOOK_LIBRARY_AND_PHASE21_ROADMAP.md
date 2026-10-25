# Authoritative Book Library Catalog (`books/`) & Phase 21+ SOTA Extraction Map

## 1. Complete Catalog of `/home/jakeshea/rust_quantlib/books/`

The workspace book library contains **86 books and monographs across 14 topic folders in `books/pdfs/`** plus **37 SOTA papers, regulatory standards, and open-access monographs in `books/deferred/21_sota_acquisitions_2026_non_china/`**.

| Folder | Representative Monographs & Authors | Core Quantitative / Institutional Domain |
|---|---|---|
| `01_quantitative_methods` | Hastie–Tibshirani–Friedman (2009) *ESL II*; Tsay (2010) *Analysis of Financial Time Series*; Tsay (2014) *Multivariate Time Series Analysis*; Greene (2019); Wooldridge (2023) | GARCH / EGARCH / GJR-GARCH, Heston–Nandi GARCH option pricing, DCC-GARCH, BEKK-GARCH, state-space Kalman filtering, Johansen cointegration |
| `02_economics` | Cochrane (2023) *Fiscal Theory of the Price Level*; Mas-Colell–Whinston–Green (1995); Romer (2019); Hart & Mas-Colell (2000) | Campbell–Shiller fiscal valuation identity, Arrow–Debreu state prices, regret-matching equilibria |
| `03_financial_reporting_analysis` | Penman (2013) *Financial Statement Analysis and Security Valuation*; Kieso et al. (2020) *IFRS*; HKFRS 9 (2024) / IAS 32 | Residual Income (`RIV`), Abnormal Earnings Growth (`AEG`), IFRS 9 Stage 1/2/3 Expected Credit Loss (`ECL`), bifurcated convertibles |
| `05_equity` | Damodaran (2025) *Investment Valuation* (4th ed.) | Multi-stage FCFF/FCFE/DDM, Implied Equity Risk Premium (`IERP`) via 1D IFT, real options |
| `06_fixed_income_and_credit` | Brigo & Mercurio (2006) *Interest Rate Models*; Crépey–Bielecki–Brigo (2014) *Counterparty Risk and Funding*; Tuckman & Serrat (2011); Veronesi (2010); Davidson & Levin (2014) *Mortgage Valuation Models*; Duffie & Singleton (2003) *Credit Risk* | LMM/SMM swaption approximations, Fong–Vasicek stochastic-vol short rate, Marshall–Olkin common-shock XVA, 4-factor mortgage prepayment & IO/PO strip OAS, Duffie–Singleton `RMV` vs. `RFV` credit basis |
| `09_portfolio_management_and_asset_pricing` | Cochrane (2005) *Asset Pricing*; Pedersen (2015) *Efficiently Inefficient* | Hansen–Jagannathan SDF bounds, Garleanu–Pedersen (2013) dynamic multi-period portfolio optimization with quadratic transaction costs & alpha decay |
| `10_behavioral_finance` | Wakker (2010) *Prospect Theory*; Barberis–Jin–Wang (2021); Shleifer (2000); Hommes & LeBaron (2018) | Cumulative Prospect Theory (`CPT`) Prelec probability weighting, DSSW noise-trader limits to arbitrage |
| `11_risk_management` | McNeil–Frey–Embrechts (2015) *Quantitative Risk Management*; Embrechts–Klüppelberg–Mikosch (1997); Föllmer & Schied (2025) *Stochastic Finance*; Potters & Bouchaud (2020) *Random Matrix Theory*; Gregory (2020) *The xVA Challenge*; Christoffersen (2012); Sornette (2017) | Peaks-Over-Threshold (`POT`) Generalized Pareto Distribution (`GPD`) extreme VaR/ES, Rotationally Invariant Estimators (`RIE`), spectral risk measures, Kupiec/Christoffersen backtesting, JLS LPPL crash hazard |
| `14_microstructure_and_trading` | Cartea–Jaimungal–Penalva (2015) *Algorithmic and High-Frequency Trading*; Foucault–Pagano–Röell (2013) *Market Liquidity*; Hasbrouck (2007); O'Hara (1995); Harris (2003) | GLFT infinite-horizon inventory market making, optimal execution with order-flow alpha, Kyle (1985) $\lambda$, Glosten–Milgrom (1985) adverse selection, EKOP `PIN`, Hasbrouck information share |
| `15_performance_and_attribution` | Bacon (2023) *Practical Portfolio Performance Measurement and Attribution*; Bacon (2022) *Risk-Adjusted Performance*; Christopherson–Carino–Ferson (2009); Colin (2016) | Brinson–Hood–Beebower (`BHB`) & Brinson–Fachler (`BF`) attribution, Carino (1999) & Menchero multi-period geometric smoothing, Karnosky–Singer multi-currency attribution, Campisi fixed-income attribution |
| `22_fund_level_arbitrage` | Madhavan (2016) *ETFs*; Hill et al. (2025); Abner (2016); van der Merwe (2015) | ETF creation/redemption AP custom-basket no-arbitrage bands, Bangia exogenous Liquidity-Adjusted VaR (`LVaR`) |
| `convertible_bonds` & `shared_anchors` | De Spiegeleer–Schoutens–Van Hulle (2014) *The Handbook of Hybrid Securities*; Calamos (2003); Koziol (2004); Thorp & Kassouf (1967); Glasserman (2003); Lando (2004); Hull (2022) | Mandatory Convertibles (`DECS` / `PERCS` / `PRIDES`), endogenous warrant dilution, make-whole & soft-call Parisian triggers, Malliavin/LRM Monte Carlo Greeks |
| `deferred/21_sota_acquisitions_2026_non_china` | Bayer–dos Reis–Horvath–Oberhauser (2025) *Signature Methods in Finance* (Springer OA); Sîrbu–Pikovsky–Shreve (2004) *Perpetual Convertible Bonds*; FEDS (2024) *Optimal Design of Contingent Capital*; López de Prado (2025) *Causality in Factor Investing* | Truncated path signatures $\mathrm{Sig}^N(X)$ & Chen's identity, structural firm-value Dynkin game convertible bonds, dual-trigger AT1 CoCos |

---

## 2. Phase 21+ SOTA Extraction Map (Zero-Overlap Institutional Candidates)

1. **Mandatory Convertibles (`DECS` / `PERCS` / `PRIDES`) & Endogenous Dilution-Adjusted Warrants** — De Spiegeleer, Schoutens & Van Hulle (2014, Ch. 11) & Koziol (2004).
2. **4-Factor Structural Mortgage Prepayment & IO/PO Strip OAS** — Davidson & Levin (2014, Ch. 6–8) & Tuckman–Serrat (2011).
3. **Fractional Recovery of Market Value (`RMV`) vs. Recovery of Face/Treasury (`RFV`/`RT`) Risky Bond & CDS-Bond Basis** — Duffie & Singleton (1999 / 2003, Ch. 5–6).
4. **Marshall–Olkin Common-Shock Bilateral CVA/DVA/FVA with Simultaneous Default Jump** — Crépey, Bielecki & Brigo (2014, Ch. 8).
5. **Path-Signature Iterated Integrals (`Chen's Identity`) & Signature Option Pricing** — Bayer, dos Reis, Horvath & Oberhauser (2025, Ch. 2–4).
6. **Asymmetric `GJR-GARCH(1,1)` / Heston–Nandi Option Pricing & `DCC-GARCH(1,1)` Dynamic Conditional Correlation** — Tsay (2010, Ch. 3; 2014, Ch. 7).
7. **Peaks-Over-Threshold (`POT`) Generalized Pareto Distribution (`GPD`) Extreme VaR/ES & Christoffersen Markov Backtesting** — McNeil, Frey & Embrechts (2015, Ch. 5) & Christoffersen (2012).
8. **Microstructure Price Impact & Informed Trading (`Kyle Lambda`, `Glosten–Milgrom`, `EKOP PIN` & `Bangia LVaR`)** — Foucault–Pagano–Röell (2013), O'Hara (1995) & van der Merwe (2015).
9. **Garleanu & Pedersen (2013) Dynamic Multi-Period Portfolio Optimization with Quadratic Transaction Costs & Mean-Reverting Alpha Decay** — Pedersen (2015, Ch. 3).
10. **Bacon (2023) Brinson–Fachler & Carino (1999) Exact Logarithmic Multi-Period Performance Attribution** — Bacon (2023, Ch. 5 & 8).

### Phase 24 Round 1 delivery
The twelve institutional engines above are now wired into `ql-fixed-income`, `ql-credit`, `ql-xva`, `ql-models`, `ql-vol`, `ql-portfolio`, and `ql-ad`, with Suite 36 (`36_mbs_io_po_rmv_mo_xva_dilution_sig_garch_pot_rie_pin_gp_brinson_ierp`) and workstation Tab 34 (`34 Phase 24 SOTA`).
