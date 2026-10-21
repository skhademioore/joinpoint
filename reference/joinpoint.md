# Fit joinpoint regression models

Fits joinpoint (segmented line) regression models to one or more trend
series (Kim et al. 2000). For each cohort (combination of `by` values)
the best fitting model with `kmin` to `kmax` joinpoints is found by an
exhaustive grid search, the number of joinpoints is selected with the
chosen criterion, and the annual percent change (APC), average annual
percent change (AAPC) and joinpoint locations are estimated with
confidence intervals.

## Usage

``` r
joinpoint(
  data,
  x,
  y = NULL,
  by = NULL,
  se = NULL,
  count = NULL,
  population = NULL,
  type = c("rate", "crude_rate", "proportion", "percent", "count", "other"),
  rate_per = 1e+05,
  family = c("gaussian", "poisson", "binomial"),
  log = TRUE,
  variance = NULL,
  autocorrelation = 0,
  vcov = NULL,
  kmin = 0,
  kmax = NULL,
  method = c("wbic", "bic", "bic3", "mbic", "dds", "permutation", "lrt"),
  ci_method = NULL,
  covariates = NULL,
  continuous = FALSE,
  overdispersion = FALSE,
  exclude = NULL,
  shift = 0,
  jump = NULL,
  comparability_ratio = NULL,
  cr_variance = 0,
  interval = "Annual",
  control = jp_control()
)
```

## Arguments

- data:

  A data frame.

- x:

  Independent variable (column name, bare or quoted), usually the
  calendar year. Values must be unique within each cohort.

- y:

  Dependent variable (a rate, proportion, count or other value). May be
  omitted when `count` and `population` are given.

- by:

  Optional character vector of by-variables defining cohorts.

- se:

  Optional column of standard errors of `y` (on the original scale).

- count, population:

  Optional columns of counts and populations (or numerators and
  denominators) used to calculate crude rates, proportions or percents
  and their standard errors.

- type:

  Type of dependent variable: `"rate"` (a provided crude or age adjusted
  rate), `"crude_rate"`, `"proportion"`, `"percent"`, `"count"` or
  `"other"`.

- rate_per:

  Rate multiplier for calculated crude rates.

- family:

  `"gaussian"` (default) fits by (weighted or generalized) least
  squares. `"poisson"` fits a Poisson regression of `count` with log
  link and `population` as exposure, and `"binomial"` a logistic
  regression of `count` events out of `population` trials, both by
  maximum likelihood (Czajkowski et al. 2008); for these the segment
  changes are percent changes in the rate or in the odds.

- log:

  Fit the log-linear model `ln(y) = xb` (`TRUE`, default) or the linear
  model `y = xb` (`FALSE`).

- variance:

  Heteroscedastic errors option: `"constant"`, `"se"` (standard errors
  provided or calculated) or `"poisson"`. Defaults to `"se"` when
  standard errors are available and `"constant"` otherwise.

- autocorrelation:

  First order autocorrelation of the errors: a number in (-1, 1), or
  `"estimate"` to estimate it from the data.

- vcov:

  Optional variance-covariance matrix of the responses (single cohort
  only).

- kmin, kmax:

  Minimum and maximum number of joinpoints. `kmax = NULL` uses the
  recommended maximum for the number of observations.

- method:

  Model selection method: `"wbic"` (weighted BIC, default), `"bic"`,
  `"bic3"`, `"mbic"` (modified BIC), `"dds"` (data dependent selection)
  or `"permutation"` for least squares fits; `"wbic"`, `"bic"`, `"bic3"`
  or `"lrt"` (sequential Monte Carlo likelihood ratio test) for Poisson
  and binomial fits.

- ci_method:

  Confidence interval method for APCs, AAPCs and joinpoint locations:
  `"empirical"` (empirical quantile, the default for least squares fits
  with grid joinpoints) or `"parametric"` (the only method for Poisson
  and binomial fits and for continuous joinpoints). Parametric intervals
  condition on the estimated joinpoints, except the interval of a
  joinpoint location; for Poisson and binomial fits they are Wald
  intervals and the joinpoint interval is based on the profile deviance.

- covariates:

  Optional one-sided formula of covariates adjusted for in every
  segment, for example `~ factor(month)` for seasonality. The covariates
  are centered, so the fitted trend refers to their average. Their names
  must differ from those of the model's own coefficients: `intercept`,
  `slope`, `jump` and `delta1`, `delta2`, and so on.

- continuous:

  Estimate joinpoints anywhere between observations (Hudson 1966)
  instead of only at the grid points: every placement of joinpoints at
  observations or inside the gaps between them is examined, so the least
  squares solution is exact. Placements are those allowed by the
  placement rules and their limits, in which an observation on a
  joinpoint counts toward one of its two segments. The limits of the
  parametric joinpoint intervals are located inside the gaps. Available
  for least squares fits with parametric intervals and selection by a
  criterion; `n_grid` is not used.

- overdispersion:

  For Poisson and binomial fits, estimate a dispersion parameter
  (quasi-likelihood) instead of fixing it at 1 (maximum likelihood): the
  Pearson statistic of the largest model divided by the number of
  observations less its coefficients and joinpoints, and at least 1. It
  scales the Wald covariances and the deviances in the criteria and is
  reported by
  [`model_selection()`](https://skhademioore.github.io/joinpoint/reference/model_selection.md).

- exclude:

  Values of `x` to exclude from model fitting. Excluded points are kept
  in the data and graphs.

- shift:

  Constant added to `x` before fitting.

- jump:

  For coding changes: the last value of `x` under the old coding. Fits
  the jump model, or the comparability ratio model when
  `comparability_ratio` is given.

- comparability_ratio, cr_variance:

  Comparability ratio (new coding over old coding) and its variance for
  the comparability ratio model.

- interval:

  Interval label of `x` used in output, e.g. `"Annual"` or `"Monthly"`.

- control:

  Settings from
  [`jp_control()`](https://skhademioore.github.io/joinpoint/reference/jp_control.md).

## Value

An object of class `joinpoint`.

## References

Kim HJ, Fay MP, Feuer EJ, Midthune DN (2000). Permutation tests for
joinpoint regression with applications to cancer rates. Statistics in
Medicine 19:335-351.

Clegg LX, Hankey BF, Tiwari R, Feuer EJ, Edwards BK (2009). Estimating
average annual percent change in trend analysis. Statistics in Medicine
28:3670-3682.

Kim HJ, Luo J, Chen HS, et al. (2017). Improved confidence interval for
average annual percent change in trend analysis. Statistics in Medicine
36:3059-3074.

Kim HJ, Chen HS, Midthune D, et al. (2023). Data-driven choice of a
model selection method in joinpoint regression. Journal of Applied
Statistics 50:1992-2013.

Hudson DJ (1966). Fitting segmented curves whose join points have to be
estimated. Journal of the American Statistical Association 61:1097-1129.

Czajkowski M, Gill R, Rempala G (2008). Model selection in logistic
joinpoint regression with applications to analyzing cohort mortality
patterns. Statistics in Medicine 27:1508-1526.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, x = year, y = rate, kmax = 3, ci_method = "parametric")
fit
#> Joinpoint regression (log-linear model, constant variance)
#> Model selection: WBIC; confidence intervals: parametric
#> 
#> Cohort: All
#>   Joinpoints: 3 (1947, 1954, 1981)
#>   1900-1947  APC = -3.720*  [-3.911, -3.528]
#>   1947-1954  APC = -14.660*  [-18.349, -10.803]
#>   1954-1981  APC = -8.967*  [-9.407, -8.526]
#>   1981-1998  APC = -4.520*  [-5.390, -3.642]
#>   AAPC Full Range (1900-1998) = -6.144*  [-6.504, -5.783]
#> 
#> * Significantly different from zero at alpha = 0.05.
apc(fit)
#>   model final  ci_method segment start  end   estimate      lower      upper
#> 1     3  TRUE parametric       1  1900 1947  -3.719816  -3.911050  -3.528200
#> 2     3  TRUE parametric       2  1947 1954 -14.659682 -18.349395 -10.803234
#> 3     3  TRUE parametric       3  1954 1981  -8.967482  -9.406536  -8.526299
#> 4     3  TRUE parametric       4  1981 1998  -4.520365  -5.390354  -3.642375
#>    statistic      p_value significant
#> 1 -37.890172 2.826091e-56        TRUE
#> 2  -7.127769 2.692260e-10        TRUE
#> 3 -38.619106 5.804778e-57        TRUE
#> 4 -10.042718 2.937842e-16        TRUE
```
