# Predictions from a joinpoint fit

Predictions include the effects of the covariates at their values in the
data (or in `newdata`). The fitted trend at the average of the
covariates, which the plots show, is the `modeled` column of
[`as.data.frame.joinpoint()`](https://skhademioore.github.io/joinpoint/reference/as.data.frame.joinpoint.md).

## Usage

``` r
# S3 method for class 'joinpoint'
predict(
  object,
  newdata = NULL,
  model = "final",
  scale = c("response", "link"),
  ...
)
```

## Arguments

- object:

  A `joinpoint` object.

- newdata:

  Optional data frame with the independent variable (and by-variables
  for multiple cohorts, and covariates of a fit with covariates).
  Defaults to the fitted data.

- model:

  `"final"` or a number of joinpoints.

- scale:

  `"response"` (original scale) or `"link"` (model scale).

- ...:

  Ignored.

## Value

A data frame with cohort keys, `x` and `fit`.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
predict(fit, newdata = data.frame(year = c(1950, 1990)))
#>      x       fit
#> 1 1950 23.159083
#> 2 1990  0.688869
```
