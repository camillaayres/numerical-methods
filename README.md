# numerical-methods
# Numerical Methods: Asian Options and Barrier Options

---

## Overview

This project applies Monte Carlo simulation, variance reduction, and closed-form sensitivity methods to two option pricing problems under geometric Brownian motion. Exercise 1 covers arithmetic and geometric Asian options; Exercise 2 covers a discretely monitored up-and-out barrier call.

All implementation is in MATLAB. The full write-up is in [`report/SMM313_Group4_report.pdf`](report/SMM313_Group4_report.pdf).

---

## Exercise 1 — Asian Options

Pricing and sensitivity analysis for three path-dependent payoffs: the arithmetic Asian call (Aₙ), its lower bound (LBₙ), and the geometric Asian call (Gₙ). Parameters: S₀ = 100, r = 0.04, σ = 0.3, T = 1, across nine (K, n) combinations (K ∈ {90, 100, 110}, n ∈ {4, 12, 50}).

### Q1.1 — Sensitivity Estimators

Likelihood ratio (LR) first-order delta estimators derived for all three payoffs. The score function is ζ(S₁|S₀)/(S₀σ√δ) where δ = T/n. A pathwise (PW) estimator for LBₙ does not exist: the indicator 1{Gₙ^(1/n) > K} introduces a jump discontinuity in S₀, violating the Lipschitz continuity required for differentiation under the expectation.

### Q1.2 — Exact Sensitivity Formulae

Closed-form delta expressions derived by differentiating the Curran (1994) formula for E(LBₙ) and the Kemna–Vorst (1990) formula for E(Gₙ) with respect to S₀. Since µₖ, σₖ, aₖ, σ̄ and T̄ do not depend on S₀, only b and d vary, giving db/dS₀ = dd/dS₀ = 1/(S₀σ̄√T̄).

### Q1.3 — Exact Numerical Results (`q3_exact_prices_sensitivities.m`)

E(LBₙ) and E(Gₙ) computed for all nine (K, n) pairs alongside their exact deltas. E(LBₙ) > E(Gₙ) throughout (AM–GM inequality), with gaps of 3–6%. Both prices decrease as K increases (further OTM) and as n increases for the parameters considered.

### Q1.4 — Monte Carlo with Control Variates (`q4_monte_carlo_control_variates.m`)

M = 10⁵ GBM paths used to estimate E(Aₙ) and dE(Aₙ)/dS₀. Two control variates compared:

- **Method (a):** LBₙ for prices, LR delta of LBₙ for sensitivities — uses Curran closed forms
- **Method (b):** Gₙ for prices, LR delta of Gₙ for sensitivities — uses Kemna–Vorst closed forms

Method (a) achieves standard errors roughly 10× smaller than method (b) for both prices and deltas. Efficiency ratios E(K,n) range from 0.009 to 0.025, meaning method (a) is 40–110× more efficient. LBₙ is a tighter control because it shares the same arithmetic path average as Aₙ, giving near-unit correlation with the target.

### Q1.5 — Convergence to Continuous Limit (`q5_convergence_continuous_limit.m`)

E(LBₙ) evaluated for n = 2², 2³, ..., 2¹⁰ and compared to the Thompson (1999) continuous-monitoring closed form E(LB∞) = 7.72696136. Successive error ratios converge to exactly 2, confirming O(1/n) convergence — each doubling of n halves the error precisely. This is faster than the O(n^{-1/2}) rate in Exercise 2, since the error here arises from quadrature approximation rather than missed Brownian bridge crossings.

---

## Exercise 2 — Discretely Monitored Up-and-Out Barrier Call

A discretely monitored up-and-out barrier call priced by Monte Carlo with antithetic variates, compared to the Shreve (2004) continuous-monitoring closed form. Parameters: S₀ = 110, K = 100, r = 0.05, σ = 0.1, T = 2, U ∈ {160, 170}, n = 2², ..., 2¹⁰, M = 10⁶ simulations per case.

### Q2.1 — Results and Accuracy (`exercise2_main.m`, `mc_uoc_antithetic.m`)

Antithetic variates reduce variance by a factor of 4–9× relative to crude MC at the same budget. The reduction is larger for U = 170 (barrier further away, payoff smoother, antithetic pairing more effective). CI widths remain roughly constant as n increases, confirming simulation error is controlled by M, not monitoring frequency. Discretisation error dominates at small n (25× the CI half-width at n = 4 for U = 160) and shrinks as n grows.

### Q2.2 — Convergence to Continuous Benchmark

Discrete prices decrease monotonically to the continuous benchmark as n increases: more monitoring dates detect more barrier crossings, knocking out more paths. Convergence follows the known O(n^{-1/2}) rate — error halves approximately every quadrupling of n. Residual gaps of 0.074 (U = 160) and 0.038 (U = 170) remain at n = 1024, consistent with the slow convergence rate. U = 170 is consistently more expensive than U = 160 (higher barrier → lower knock-out probability).

---

## Key Results

| | E(LBₙ) | E(Gₙ) | E(Aₙ) MC (method a) |
|---|---|---|---|
| K=100, n=12 | 8.2345 | 7.8021 | 8.2382 ± 0.000166 |
| K=100, n=50 | 7.8490 | 7.4155 | 7.8523 ± 0.000152 |

| Barrier | n | MC Price | Gap to CF | Var. Ratio (AV vs crude) |
|---|---|---|---|---|
| U=160 | 4 | 18.6133 | 0.659 | 4.65 |
| U=160 | 1024 | 18.0279 | 0.074 | — |
| U=170 | 4 | 19.5144 | 0.275 | 8.94 |
| U=170 | 1024 | 19.2772 | 0.038 | — |

---

## Repository Structure

```
├── README.md
├── report/
│   └── SMM313_Group4_report.pdf
├── code/
│   ├── exercise1/
│   │   ├── q3_exact_prices_sensitivities.m
│   │   ├── q4_monte_carlo_control_variates.m
│   │   ├── q5_convergence_continuous_limit.m
│   │   └── fig_pw_discontinuity.m
│   └── exercise2/
│       ├── exercise2_main.m
│       ├── mc_uoc_antithetic.m
│       ├── mc_uoc_crude.m
│       └── fUOC_continuous.m
└── results/
    └── figures/
        ├── q1_pw_discontinuity_lbn.pdf
        ├── q5_lbn_convergence_continuous_limit.pdf
        ├── q2_uoc_convergence_to_continuous.pdf
        ├── q2_uoc_discretisation_gap.pdf
        └── q2_uoc_ci_width_vs_monitoring_dates.pdf
```

---

## Requirements

MATLAB (tested on R2024b). No additional toolboxes required beyond the base installation. All scripts are self-contained.

---

## Key References

- Curran, M. (1994). Valuing Asian and portfolio options by conditioning on the geometric mean price. *Management Science*, 40.
- Kemna, A.G.Z. and Vorst, A.C.F. (1990). A pricing method for options based on average asset values. *Journal of Banking & Finance*, 14.
- Shreve, S. (2004). *Stochastic Calculus for Finance II: Continuous-Time Models*. Springer Finance.
- Thompson, G.W.P. (1999). Fast narrow bounds on the value of Asian options. Working paper, Cambridge University.
