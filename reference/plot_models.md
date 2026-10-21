# Plot all fitted models

The best fitting model for every number of joinpoints of one cohort,
either as small multiples (the selection statistic in each panel title)
or overlaid in one panel, with the final selected model highlighted.

## Usage

``` r
plot_models(object, cohort = 1L, overlay = FALSE, ncol = NULL)
```

## Arguments

- object:

  A `joinpoint` object.

- cohort:

  Index of the cohort to draw.

- overlay:

  Draw all models in one panel.

- ncol:

  Number of columns of the small multiples.

## Value

A ggplot object.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
plot_models(fit)
```
