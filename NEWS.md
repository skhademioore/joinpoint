# joinpoint 0.1.0

## Models

* `joinpoint()` fits joinpoint regression models to one or more series
  (cohorts given by `by`): linear or log-linear models of rates, age adjusted
  rates, proportions, percents, counts or other values, with constant
  variance, supplied standard errors, Poisson variance, first order
  autocorrelated errors (estimated at the actual gaps when observations are
  excluded) or a supplied covariance matrix.
* Joinpoints are found by an exhaustive grid search with placement
  constraints and optional points between observations (Lerman 1980),
  compiled in C++ and evaluated from cross products, with an exact search
  from residual vectors when a model fits almost perfectly. Fitted values
  are evaluated in centered coordinates, so shifted or distant x values keep
  full precision.
* Poisson and logistic joinpoint models fitted by exact maximum likelihood
  (Czajkowski et al. 2008), with every candidate compared at convergence,
  Wald intervals from the QR factorization, profile deviance joinpoint
  intervals, deviance based criteria and optional quasi-likelihood dispersion.
* Continuous joinpoints between observations by an exact search (Hudson 1966).
* Covariates (`covariates = ~ ...`), excluded data points, shifted x values,
  jump and comparability ratio models for coding changes (Chen et al. 2020),
  with the uncertainty of the comparability ratio shared by all old coding
  observations (its known variance kept apart from the estimated residual
  variance) or added to each of them as in the NCI Joinpoint program.

## Model selection

* WBIC, BIC, BIC3, modified BIC, data dependent selection, the sequential
  permutation test (Kim et al. 2000) and a sequential Monte Carlo likelihood
  ratio test, with optional early stopping (Besag and Clifford 1991) and a
  minimum APC difference worth detecting.
* A time budget lowers the recommended maximum number of joinpoints when the
  exhaustive search would be too slow, also across the age group blocks of
  clustering.

## Inference

* APC, AAPC (full range, last N observations, custom ranges) and joinpoint
  location estimates with parametric or empirical quantile confidence
  intervals and p values (Kim et al. 2017).
* Pairwise tests of parallelism and coincidence with AAPC comparisons and
  multiplicity adjustment (`jp_compare()`; Kim et al. 2004).
* Clustering of contiguous age groups (`jp_cluster()`), with final cluster
  rates standardized by the original or the compressed age groups (the
  settings are kept in the result), and
  direct age adjustment with delay factors (`age_adjust()`).

## Output and graphics

* Exact fits are recognized: a constant series has a zero trend with p = 1,
  an unidentified joinpoint gets the whole admissible interval, and
  parameters that are not estimable get no standard error.
* Accessors `apc()`, `aapc()`, `joinpoints()`, `model_selection()`,
  `model_estimates()`, `jump_estimates()`, the publication table `jp_table()`,
  `coef()` and `vcov()` methods with matching parameters, `fitted()`,
  `residuals()`, `predict()` (one prediction per row of `newdata`) and
  `as.data.frame()` methods, and `jp_export()` for delimited text files.
* `ggplot2` graphics: `autoplot()` trend figures with segment labels placed to
  avoid each other, the curves, the points and their error bars (exact
  Poisson or binomial limits for rates calculated from counts);
  `plot_apc()`, `plot_aapc()`, `plot_models()`, `plot_selection()`,
  `plot_residuals()` (autocorrelation lags counted in time steps across
  excluded observations), comparison plots, `theme_joinpoint()` and color
  scales.
* Interactive graphics (`jp_interactive()`, or `interactive = TRUE` in the
  plot functions) built with ggiraph: tooltips with the APC, confidence
  interval and p value of each segment, the values of the observations and
  the joinpoint intervals, linked highlighting and zoom, for the RStudio
  viewer, HTML documents, pkgdown sites and Shiny.

## Export and reports

* `jp_save()` writes tables to delimited text, Markdown, LaTeX (booktabs),
  Word, OpenDocument, PowerPoint, RTF, HTML and PDF; figures to the image
  formats of ggplot2, LaTeX figures and documents; and interactive graphics
  to standalone web pages.
* `jp_report()` writes a complete report of a fit as an R Markdown or Quarto
  source, or renders it to HTML, Word, OpenDocument, PDF, Markdown or LaTeX.
* With `options(joinpoint.knit_table = TRUE)`, fits printed in R Markdown
  and Quarto documents appear as formatted tables.

## Data and documentation

* Data sets `us_mortality`, `std_populations` and `ontario_bmi`.
* Nine vignettes, including a reproduction of a published analysis of
  prepregnancy body mass index (Khademioore et al. 2026).
