# Observed and modeled values

Observed values, standard errors, the supplied counts and populations,
modeled values (the fitted trend, with covariates at their average) for
the selected model (or all models) and flags for excluded points, zero
counts replaced by 0.5 in the observed value and joinpoints.

## Usage

``` r
# S3 method for class 'joinpoint'
as.data.frame(x, row.names = NULL, optional = FALSE, model = "final", ...)
```

## Arguments

- x:

  A `joinpoint` object.

- row.names, optional:

  Ignored.

- model:

  `"final"`, a number of joinpoints, or `"all"`.

- ...:

  Ignored.

## Value

A data frame.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
head(as.data.frame(fit))
#>   model final    x observed se count population  modeled included zero_count
#> 1     2  TRUE 1900    206.0 NA    NA         NA 237.7315     TRUE      FALSE
#> 2     2  TRUE 1901    200.9 NA    NA         NA 228.9907     TRUE      FALSE
#> 3     2  TRUE 1902    183.8 NA    NA         NA 220.5712     TRUE      FALSE
#> 4     2  TRUE 1903    185.4 NA    NA         NA 212.4613     TRUE      FALSE
#> 5     2  TRUE 1904    195.4 NA    NA         NA 204.6496     TRUE      FALSE
#> 6     2  TRUE 1905    187.3 NA    NA         NA 197.1251     TRUE      FALSE
#>   joinpoint segment       apc
#> 1     FALSE       1 -3.676767
#> 2     FALSE       1 -3.676767
#> 3     FALSE       1 -3.676767
#> 4     FALSE       1 -3.676767
#> 5     FALSE       1 -3.676767
#> 6     FALSE       1 -3.676767
```
