# Jump and comparability ratio estimates

For jump models (see the `jump` argument of
[`joinpoint()`](https://skhademioore.github.io/joinpoint/reference/joinpoint.md)),
the estimated jump on the log scale and the comparability ratio
`exp(jump)` with standard errors, confidence intervals and tests (null
values 0 and 1).

## Usage

``` r
jump_estimates(object, model = "final")
```

## Arguments

- object:

  A `joinpoint` object.

- model:

  `"final"` for the final selected model, a number of joinpoints, or
  `"all"` for every fitted model.

## Value

A data frame, or `NULL` when no jump model was fitted.

## Examples

``` r
heart <- subset(us_mortality, cause == "Heart Disease" & year >= 1950)
fit <- joinpoint(heart, year, rate, kmax = 2, jump = 1978,
                 ci_method = "parametric")
jump_estimates(fit)
#>   model final  ci_method           statistic   estimate         se       lower
#> 1     1  TRUE parametric          Jump value 0.03430509 0.01272229 0.008648136
#> 2     1  TRUE parametric Comparability ratio 1.03490029 0.01316630 1.009094808
#>        upper test_statistic     p_value
#> 1 0.05996203       2.696455 0.009965438
#> 2 1.06070577       2.650728 0.008031852
```
