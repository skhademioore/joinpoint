# Publication table of trends

Formats the APCs of the segments (slopes for linear models) and the
AAPCs (AACs) of the final, or a given, model as a character table with
estimates, confidence intervals and p values, ready for
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) or other
table packages.

## Usage

``` r
jp_table(object, model = "final", digits = 2, aapc = TRUE)
```

## Arguments

- object:

  A `joinpoint` object.

- model:

  `"final"` for the final selected model, a number of joinpoints, or
  `"all"` for every fitted model.

- digits:

  Decimals of estimates and confidence limits.

- aapc:

  Append the AAPC rows.

## Value

A data frame of character columns.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 3, ci_method = "parametric")
jp_table(fit)
#>       Trend    Period      APC / AAPC (95% CI)       p
#> 1 Segment 1 1900-1947    -3.72* (-3.91, -3.53) < 0.001
#> 2 Segment 2 1947-1954 -14.66* (-18.35, -10.80) < 0.001
#> 3 Segment 3 1954-1981    -8.97* (-9.41, -8.53) < 0.001
#> 4 Segment 4 1981-1998    -4.52* (-5.39, -3.64) < 0.001
#> 5      AAPC 1900-1998    -6.14* (-6.50, -5.78) < 0.001
```
