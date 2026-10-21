# Pairwise comparison of joinpoint trends

Tests whether the joinpoint regression functions of two cohorts are
identical (test of coincidence) or parallel (test of parallelism, equal
joinpoints and slopes with different intercepts), using the permutation
test of Kim et al. (2004). Every pair of levels of the innermost
by-variable is compared within each combination of the other
by-variables. The test is carried out with `K = max(K0, K1, K2)`
joinpoints, where `K1` and `K2` are the numbers selected for each cohort
and `K0` the number selected for the two cohorts fitted together under
the null hypothesis. For each pair whose test rejects the null
hypothesis (after any p value adjustment) the AAPCs of the two cohorts
of a log-linear fit are compared by a z test with parametric variances,
also when the fit used empirical quantile intervals; these follow-up p
values are not adjusted.

## Usage

``` r
jp_compare(
  object,
  type = c("parallel", "coincidence"),
  alpha = 0.05,
  n_perm = 4499,
  seed = 7160,
  p_adjust = "none"
)
```

## Arguments

- object:

  A `joinpoint` object with at least one by-variable, fitted with
  uncorrelated errors.

- type:

  `"parallel"` or `"coincidence"`.

- alpha:

  Significance level of the test.

- n_perm:

  Number of permuted data sets.

- seed:

  Random number seed.

- p_adjust:

  Method of [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)
  applied to the p values of the tests of all pairs, for example
  `"holm"` when many pairs are compared; the unadjusted p values are
  kept in `p_unadjusted`.

## Value

An object of class `jp_comparison`: a data frame with one row per pair
and the attributes `aapc` (AAPC comparisons), `fits` (null model fits)
and `p_adjust`.

## References

Kim HJ, Fay MP, Yu B, Barrett MJ, Feuer EJ (2004). Comparability of
segmented line regression models. Biometrics 60:1005-1014.

## Examples

``` r
d <- subset(us_mortality, cause %in% c("Stroke", "Heart Disease"))
fit <- joinpoint(d, year, rate, by = "cause", kmax = 2,
                 ci_method = "parametric")
cmp <- jp_compare(fit, n_perm = 199)
cmp
#> Pairwise comparison: test of parallel
#> 
#>         group1 group2     test joinpoints k_null k_group1 k_group2 statistic
#>  Heart Disease Stroke parallel          2      1        2        2     732.1
#>  df1 df2 p_value alpha       result
#>    5 186   0.005  0.05 Not parallel
#> 
#> AAPC comparisons (group 1 minus group 2):
#>         group1 group2                  range    aapc1  aapc2 difference lower
#>  Heart Disease Stroke Full Range (1900-1998) -0.02318 -1.472      1.449 1.255
#>  upper statistic p_value significant  ci_method
#>  1.643     14.55 5.5e-48        TRUE parametric
```
