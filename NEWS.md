# joinpoint 0.1.0

First version.

## Models

* `joinpoint()` fits joinpoint regression models to one or more series
  (cohorts given by `by`): linear or log-linear models of rates, age adjusted
  rates, proportions, percents, counts or other values, with constant
  variance, supplied standard errors, Poisson variance, first order
  autocorrelated errors or a supplied covariance matrix.
* Joinpoints are found by an exhaustive grid search with placement
  constraints and optional points between observations (Lerman 1980),
  compiled in C++ and evaluated from cross products.
* Poisson and logistic joinpoint models fitted by exact maximum likelihood
  (Czajkowski et al. 2008), with Wald intervals, profile deviance joinpoint
  intervals, deviance based criteria and optional quasi-likelihood dispersion.
* Continuous joinpoints between observations by an exact search (Hudson 1966).
* Covariates (`covariates = ~ ...`), excluded data points, shifted x values,
  jump and comparability ratio models for coding changes (Chen et al. 2020).

## Model selection

* WBIC, BIC, BIC3, modified BIC, data dependent selection, the sequential
  permutation test (Kim et al. 2000) and a sequential Monte Carlo likelihood
  ratio test, with optional early stopping (Besag and Clifford 1991) and a
  minimum APC difference worth detecting.
* A time budget lowers the recommended maximum number of joinpoints when the
  exhaustive search would be too slow.

## Inference

* APC, AAPC (full range, last N observations, custom ranges) and joinpoint
  location estimates with parametric or empirical quantile confidence
  intervals and p values (Kim et al. 2017).
* Pairwise tests of parallelism and coincidence with AAPC comparisons and
  multiplicity adjustment (`jp_compare()`; Kim et al. 2004).
* Clustering of contiguous age groups (`jp_cluster()`) and direct age
  adjustment with delay factors (`age_adjust()`).

## Output and graphics

* Accessors `apc()`, `aapc()`, `joinpoints()`, `model_selection()`,
  `model_estimates()`, `jump_estimates()`, the publication table `jp_table()`,
  `coef()`, `vcov()`, `fitted()`, `residuals()`, `predict()` and
  `as.data.frame()` methods, and `jp_export()` for delimited text files.
* `ggplot2` graphics: `autoplot()` trend figures with segment labels placed to
  avoid each other and the curves, `plot_apc()`, `plot_aapc()`,
  `plot_models()`, `plot_selection()`, `plot_residuals()`, comparison plots,
  `theme_joinpoint()` and color scales.

## Data and documentation

* Data sets `us_mortality`, `std_populations` and `ontario_bmi`.
* Eight vignettes, including a reproduction of a published analysis of
  prepregnancy body mass index (Khademioore et al. 2026).
