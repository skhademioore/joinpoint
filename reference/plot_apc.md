# Plot segment APCs and AAPCs

`plot_apc()` draws, for each cohort, the APC of every segment as a step
over the years it covers, with its confidence interval as a shaded band
and a zero reference line, so changes in trend read directly off the
time axis. `plot_aapc()` draws a forest plot of AAPCs with confidence
intervals by cohort and range.

## Usage

``` r
plot_apc(object, model = "final", ncol = NULL)

plot_aapc(object, model = "final")
```

## Arguments

- object:

  A `joinpoint` object.

- model:

  `"final"` or a number of joinpoints.

- ncol:

  Number of facet columns.

## Value

A ggplot object.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
plot_apc(fit)

plot_aapc(fit)
```
