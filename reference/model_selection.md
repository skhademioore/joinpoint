# Model selection results

Returns, for each cohort, the selection statistics of all fitted models
(SSE, BIC, BIC3 and, depending on the method, the WBIC weight and value
or the modified BIC; for Poisson and binomial fits the deviance and
dispersion, the criteria using the deviance divided by the dispersion),
R squared and the Durbin-Watson statistic (undefined for constant series
and exact fits), the permutation or likelihood ratio tests performed,
the data dependent selection details and the MADWD steps.

## Usage

``` r
model_selection(object)
```

## Arguments

- object:

  A `joinpoint` object.

## Value

A list with data frames `criteria`, `tests`, `dds` and `madwd` (the last
three may be `NULL`), and `selected` with the selected number of
joinpoints per cohort.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
model_selection(fit)$criteria
#>   joinpoints       sse       bic      bic3    weight      wbic df        mse
#> 1          0 14.763921 -1.810103 -1.810103        NA -1.810103 97 0.15220537
#> 2          1  3.220021 -3.240071 -3.193655 0.7818994 -3.203779 95 0.03389496
#> 3          2  1.068503 -4.250369 -4.157538 0.9336905 -4.163694 93 0.01148928
#>   r_squared durbin_watson selected
#> 1 0.9645795    0.03403132    FALSE
#> 2 0.9922748    0.13868568    FALSE
#> 3 0.9974365    0.37665278     TRUE
```
