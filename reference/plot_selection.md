# Plot model selection

Selection statistic against the number of joinpoints for each cohort.
For the permutation test the sequence of tests is drawn instead: the p
value of each test against its significance level.

## Usage

``` r
plot_selection(object, interactive = FALSE)
```

## Arguments

- object:

  A `joinpoint` object.

- interactive:

  Draw interactive `ggiraph` layers (tooltips and linked highlighting);
  show the result with
  [`jp_interactive()`](https://skhademioore.github.io/joinpoint/reference/jp_interactive.md).

## Value

A ggplot object.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
plot_selection(fit)
```
