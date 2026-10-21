# Residual diagnostics

Standardized residuals of the selected model (weighted by the error
structure, so constant variance is expected) against the independent
variable, as a normal quantile plot, or as an autocorrelation function
to judge whether an autocorrelated errors model is warranted.

## Usage

``` r
plot_residuals(object, type = c("time", "qq", "acf"), model = "final")
```

## Arguments

- object:

  A `joinpoint` object.

- type:

  `"time"`, `"qq"` or `"acf"`.

- model:

  `"final"` or a number of joinpoints.

## Value

A ggplot object.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
plot_residuals(fit)

plot_residuals(fit, type = "acf")
```
