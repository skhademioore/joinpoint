# Plot a joinpoint fit

Publication quality graph of observed values and the fitted joinpoint
model: observed points (excluded points hollow), the fitted curve,
joinpoints, and the APC of each segment as a direct label (`*` marks
APCs significantly different from zero). A coding change is marked by a
dotted line with its comparability ratio.

## Usage

``` r
# S3 method for class 'joinpoint'
autoplot(
  object,
  model = "final",
  color_by = c("cohort", "trend"),
  layout = c("facet", "overlay"),
  cohorts = NULL,
  labels = NULL,
  aapc = FALSE,
  band = TRUE,
  joinpoint_ci = FALSE,
  error_bars = FALSE,
  log_y = FALSE,
  digits = NULL,
  ncol = NULL,
  ...
)

# S3 method for class 'joinpoint'
plot(x, ...)
```

## Arguments

- object:

  A `joinpoint` object.

- model:

  `"final"` or a number of joinpoints.

- color_by:

  `"cohort"` colors each cohort; `"trend"` colors each segment by its
  direction: significantly increasing, significantly decreasing, stable,
  or not tested.

- layout:

  `"facet"` draws one panel per cohort; `"overlay"` draws all cohorts in
  one panel.

- cohorts:

  Optional integer indices of cohorts to draw.

- labels:

  Segment labels: `"apc"` (APC, or slope for linear models), `"ci"` (APC
  with its confidence interval) or `"none"`. The default is `"apc"`, or
  `"none"` when several cohorts are overlaid.

- aapc:

  Add the AAPC of the first range (the full range by default, see
  `aapc_full` in
  [`jp_control()`](https://skhademioore.github.io/joinpoint/reference/jp_control.md))
  to the title (one cohort), to each panel title (facets) or to each
  legend entry (overlaid cohorts).

- band:

  Draw the pointwise 95% band of the fitted curve from the empirical
  quantile resamples, when available.

- joinpoint_ci:

  Shade the confidence interval of each joinpoint location.

- error_bars:

  Draw 95% confidence intervals of the observed values (observed value
  plus or minus 1.96 standard errors) when standard errors are
  available.

- log_y:

  Use a logarithmic y axis.

- digits:

  Decimals of the segment labels (default 1 for APCs and 3 for slopes of
  linear models).

- ncol:

  Number of facet columns.

- ...:

  Ignored.

- x:

  A `joinpoint` object.

## Value

A ggplot object.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 3, ci_method = "parametric")
ggplot2::autoplot(fit)

ggplot2::autoplot(fit, color_by = "trend", log_y = TRUE)
```
