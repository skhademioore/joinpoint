# Model estimates

Regression coefficients of the selected (or a specified) model in the
general changepoint parameterization: segment intercepts, segment
slopes, slope changes and any jump and covariate coefficients, with
standard errors from the unconstrained fit (Wald standard errors for
Poisson and binomial fits), test statistics and p values.

## Usage

``` r
model_estimates(object, model = "final")
```

## Arguments

- object:

  A `joinpoint` object.

- model:

  `"final"` for the final selected model, a number of joinpoints, or
  `"all"` for every fitted model.

## Value

A data frame.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
model_estimates(fit)
#>   model final  ci_method         parameter     estimate          se  statistic
#> 1     2  TRUE parametric       Intercept 1  76.64635549 2.550114550  30.056044
#> 2     2  TRUE parametric       Intercept 2 203.14199030 3.399513282  59.756198
#> 3     2  TRUE parametric       Intercept 3  86.82396790 9.763743918   8.892487
#> 4     2  TRUE parametric           Slope 1  -0.03746064 0.001327466 -28.219669
#> 5     2  TRUE parametric           Slope 2  -0.10256390 0.001733095 -59.179622
#> 6     2  TRUE parametric           Slope 3  -0.04381742 0.004907620  -8.928446
#> 7     2  TRUE parametric Slope 2 - Slope 1  -0.06510326 0.002183067 -29.821921
#> 8     2  TRUE parametric Slope 3 - Slope 2   0.05874648 0.005204647  11.287311
#>        p_value df
#> 1 4.893626e-49 91
#> 2 8.219255e-75 91
#> 3 5.301410e-14 91
#> 4 8.755056e-47 91
#> 5 1.943137e-74 91
#> 6 4.458435e-14 91
#> 7 9.343411e-49 91
#> 8 5.433400e-19 91
```
