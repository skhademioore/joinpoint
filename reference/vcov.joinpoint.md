# Covariance of the segment parameters

Covariance matrix of the segment intercepts and slopes and of any jump
and covariate coefficients, from the unconstrained segment fit used for
the standard errors (the Wald covariance for Poisson and binomial fits).

## Usage

``` r
# S3 method for class 'joinpoint'
vcov(object, model = "final", ...)
```

## Arguments

- object:

  A `joinpoint` object.

- model:

  `"final"` or a number of joinpoints.

- ...:

  Ignored.

## Value

A matrix, or a list of matrices named by cohort.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
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
```
