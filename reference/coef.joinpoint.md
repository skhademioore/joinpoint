# Coefficients and their covariance

[`coef()`](https://rdrr.io/r/stats/coef.html) returns, by default, the
segment intercepts and slopes followed by any jump and covariate
coefficients, the parameters whose covariance
[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns: names and order
agree. The covariance is that of the unconstrained segment fit used for
the standard errors (the Wald covariance for Poisson and binomial fits).
`type = "hinge"` gives instead the coefficients of the continuous fit in
its hinge form: intercept, first slope and the slope changes `delta1`,
`delta2`, ... at the joinpoints.

## Usage

``` r
# S3 method for class 'joinpoint'
coef(object, model = "final", type = c("segment", "hinge"), ...)

# S3 method for class 'joinpoint'
vcov(object, model = "final", ...)
```

## Arguments

- object:

  A `joinpoint` object.

- model:

  `"final"` or a number of joinpoints.

- type:

  `"segment"` or `"hinge"` (for
  [`coef()`](https://rdrr.io/r/stats/coef.html)).

- ...:

  Ignored.

## Value

A vector or matrix, or a list of them named by cohort.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
coef(fit)
#>  Intercept 1      Slope 1  Intercept 2      Slope 2  Intercept 3      Slope 3 
#>  76.64635549  -0.03746064 203.14199030  -0.10256390  86.82396790  -0.04381742 
vcov(fit)
#>              Intercept 1       Slope 1  Intercept 2       Slope 2 Intercept 3
#> Intercept 1  6.503084219 -3.385119e-03  0.000000000  0.000000e+00  0.00000000
#> Slope 1     -0.003385119  1.762165e-06  0.000000000  0.000000e+00  0.00000000
#> Intercept 2  0.000000000  0.000000e+00 11.556690555 -5.891596e-03  0.00000000
#> Slope 2      0.000000000  0.000000e+00 -0.005891596  3.003618e-06  0.00000000
#> Intercept 3  0.000000000  0.000000e+00  0.000000000  0.000000e+00 95.33069530
#> Slope 3      0.000000000  0.000000e+00  0.000000000  0.000000e+00 -0.04791659
#>                   Slope 3
#> Intercept 1  0.000000e+00
#> Slope 1      0.000000e+00
#> Intercept 2  0.000000e+00
#> Slope 2      0.000000e+00
#> Intercept 3 -4.791659e-02
#> Slope 3      2.408474e-05
coef(fit, type = "hinge")
#>   intercept       slope      delta1      delta2 
#> 76.64635549 -0.03746064 -0.06510326  0.05874648 
```
