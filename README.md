------------------------------------------------------------------------

editor_options: markdown: wrap: 72 ---

# MDPL

This package is the software implementation of the MDPL research paper (Multi-Dimensional Power-Likelihood), a quasi-empirical Bayesian framework that jointly calibrates borrowing from multiple auxiliary datasets when estimating treatment effects in a randomised trial.

The likelihood of each auxiliary (e.g. observational) data set is raised to a weight chosen from the trial data, so informative data sets are borrowed from and biased ones are down-weighted. The package returns conditional (CATE) and average (ATE) treatment effects, with optional bootstrap intervals that account for the weight selection.

## Quick start

``` r
library(MDPL)

# the trial, then a named list of observational data sets
fit <- mdpl(Y ~ A * X, trial, list(src1, src2), treat = "A")
fit                                       # weights, coefficients, ATE and CATE formula
fit$eta                                   # one borrowing weight per source

fit$estimates$ate                         # average treatment effect (trial covariates)
fit$estimates$cate(data.frame(X = 0:2))   # CATE, predicted on new covariates

mdpl_cate(fit, data.frame(X = 0:2))       # recompute for other covariates or scales
mdpl_ate(fit)                             # ... or populations, via `data =`

# bootstrap the whole procedure (refit and reselect weights): off by default
fit <- mdpl(Y ~ A * X, trial, list(src1, src2), treat = "A",
            bootstrap_CI = TRUE, B = 500)
```

The fit stores the weights (`fit$eta`) and, in `fit$estimates`, the treatment-effect coefficients (`coef`), the ATE over the trial covariates (`ate`) and the CATE as a formula and function (`cate`). Each can be switched off in `mdpl()`.

## Models

| specification | model |
|----|----|
| `mdpl_glm(formula, family)` (or just a formula) | Gaussian, binomial, Poisson, Gamma GLMs |
| `mdpl_causl(formulas, family)` | frugal Gaussian-copula models from the `causl` package (suggested) |
| `mdpl_custom(loglik, start, ...)` | any smooth parametric likelihood you supply; check it with `mdpl_check_model()` |

Models must be regular parametric likelihoods with the same parameterisation in the trial and every source; trees, penalised fits and nonparametric models are out of scope, and every source must target the same quantity as the trial. Survival outcomes are in the companion package `gbi`. See the `mdpl-models` vignette.

## Weight selection

Each source gets one weight, chosen from the trial data and starting from no borrowing.

- `method = "taylor"` (default): fast, maximises a quadratic surrogate of the trial's predictive score; handles any number of sources.
- `method = "elpd"`: an ELPD proxy (WAIC of the trial under the power posterior).
- `method = "cv"`: cross-validation of the fused estimate on held-out trial units; slowest.

`target = "cate"` judges the weights on the treatment-effect terms only, and `pooled = TRUE` gives all sources one common weight. See the `mdpl-weights` vignette.

## Vignettes

- `mdpl`: a first analysis, end to end.
- `mdpl-models`: GLM, causl and custom-likelihood specifications.
- `mdpl-weights`: how the borrowing weights are chosen, and how they behave.
- `mdpl-intervals`: posterior and bootstrap intervals, and when to use which.
- `mdpl-overlap`: poor covariate overlap, parametric extrapolation, and why the target stays the trial population.
- `mdpl-frugal-cate`: CATE with the frugal (copula) model and a hidden confounder.

## Status

Version 0.1.0, first cut. Planned: ordinal and multinomial outcomes, clustered data, a simultaneous band for the CATE curve, and MCMC for the power posterior.
