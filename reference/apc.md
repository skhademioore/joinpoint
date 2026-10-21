# Trend estimates from a joinpoint fit

`apc()` returns the annual percent change of each segment (the slope of
each segment for linear models), `aapc()` the average annual percent
change over the requested ranges, and `joinpoints()` the estimated
joinpoint locations with confidence intervals.

## Usage

``` r
apc(object, model = "final")

aapc(object, model = "final", ranges = NULL)

joinpoints(object, model = "final")
```

## Arguments

- object:

  A `joinpoint` object.

- model:

  `"final"` for the final selected model, a number of joinpoints, or
  `"all"` for every fitted model.

- ranges:

  Optional list of `c(start, end)` ranges; recomputes the AAPC with
  parametric intervals for these ranges, clipped to the included data.

## Value

A data frame with one row per segment, range or joinpoint and the cohort
keys.

## Details

For linear models (`log = FALSE`) `aapc()` returns the average annual
change (AAC), the length weighted mean of the segment slopes over the
range, with the same interval construction as the AAPC.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 3, ci_method = "parametric")
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
aapc(fit)
#>   model final  ci_method                  range start  end  estimate
#> 1     3  TRUE parametric Full Range (1900-1998)  1900 1998 -6.144387
#>   within_segment    lower     upper statistic       p_value      se_log
#> 1          FALSE -6.50424 -5.783149 -32.35379 1.227243e-229 0.001959975
#>   log_estimate significant
#> 1  -0.06341262        TRUE
joinpoints(fit)
#>   model final  ci_method joinpoint estimate lower upper
#> 1     3  TRUE parametric         1     1947  1913  1948
#> 2     3  TRUE parametric         2     1954  1944  1957
#> 3     3  TRUE parametric         3     1981  1979  1984
```
