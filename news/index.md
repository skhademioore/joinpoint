# Changelog

## joinpoint 0.1.0

First version.

### Models

- [`joinpoint()`](https://skhademioore.github.io/joinpoint/reference/joinpoint.md)
  fits joinpoint regression models to one or more series (cohorts given
  by `by`): linear or log-linear models of rates, age adjusted rates,
  proportions, percents, counts or other values, with constant variance,
  supplied standard errors, Poisson variance, first order autocorrelated
  errors or a supplied covariance matrix.
- Joinpoints are found by an exhaustive grid search with placement
  constraints and optional points between observations (Lerman 1980),
  compiled in C++ and evaluated from cross products.
- Poisson and logistic joinpoint models fitted by exact maximum
  likelihood (Czajkowski et al. 2008), with Wald intervals, profile
  deviance joinpoint intervals, deviance based criteria and optional
  quasi-likelihood dispersion.
- Continuous joinpoints between observations by an exact search (Hudson
  1966).
- Covariates (`covariates = ~ ...`), excluded data points, shifted x
  values, jump and comparability ratio models for coding changes (Chen
  et al. 2020).

### Model selection

- WBIC, BIC, BIC3, modified BIC, data dependent selection, the
  sequential permutation test (Kim et al. 2000) and a sequential Monte
  Carlo likelihood ratio test, with optional early stopping (Besag and
  Clifford 1991) and a minimum APC difference worth detecting.
- A time budget lowers the recommended maximum number of joinpoints when
  the exhaustive search would be too slow.

### Inference

- APC, AAPC (full range, last N observations, custom ranges) and
  joinpoint location estimates with parametric or empirical quantile
  confidence intervals and p values (Kim et al. 2017).
- Pairwise tests of parallelism and coincidence with AAPC comparisons
  and multiplicity adjustment
  ([`jp_compare()`](https://skhademioore.github.io/joinpoint/reference/jp_compare.md);
  Kim et al. 2004).
- Clustering of contiguous age groups
  ([`jp_cluster()`](https://skhademioore.github.io/joinpoint/reference/jp_cluster.md))
  and direct age adjustment with delay factors
  ([`age_adjust()`](https://skhademioore.github.io/joinpoint/reference/age_adjust.md)).

### Output and graphics

- Accessors
  [`apc()`](https://skhademioore.github.io/joinpoint/reference/apc.md),
  [`aapc()`](https://skhademioore.github.io/joinpoint/reference/apc.md),
  [`joinpoints()`](https://skhademioore.github.io/joinpoint/reference/apc.md),
  [`model_selection()`](https://skhademioore.github.io/joinpoint/reference/model_selection.md),
  [`model_estimates()`](https://skhademioore.github.io/joinpoint/reference/model_estimates.md),
  [`jump_estimates()`](https://skhademioore.github.io/joinpoint/reference/jump_estimates.md),
  the publication table
  [`jp_table()`](https://skhademioore.github.io/joinpoint/reference/jp_table.md),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  methods, and
  [`jp_export()`](https://skhademioore.github.io/joinpoint/reference/jp_export.md)
  for delimited text files.
- `ggplot2` graphics:
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  trend figures with segment labels placed to avoid each other and the
  curves,
  [`plot_apc()`](https://skhademioore.github.io/joinpoint/reference/plot_apc.md),
  [`plot_aapc()`](https://skhademioore.github.io/joinpoint/reference/plot_apc.md),
  [`plot_models()`](https://skhademioore.github.io/joinpoint/reference/plot_models.md),
  [`plot_selection()`](https://skhademioore.github.io/joinpoint/reference/plot_selection.md),
  [`plot_residuals()`](https://skhademioore.github.io/joinpoint/reference/plot_residuals.md),
  comparison plots,
  [`theme_joinpoint()`](https://skhademioore.github.io/joinpoint/reference/theme_joinpoint.md)
  and color scales.

### Data and documentation

- Data sets `us_mortality`, `std_populations` and `ontario_bmi`.
- Eight vignettes, including a reproduction of a published analysis of
  prepregnancy body mass index (Khademioore et al. 2026).
