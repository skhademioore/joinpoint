# Directly age adjusted rates

Computes directly age standardized rates and their standard errors from
age specific counts and populations. With delay factors the result is
the delay adjusted age adjusted rate.

## Usage

``` r
age_adjust(
  data,
  x,
  count,
  population,
  age,
  standard,
  by = NULL,
  delay = NULL,
  rate_per = 1e+05
)
```

## Arguments

- data:

  A data frame with one row per age group, time point and cohort.

- x:

  Column of the independent variable (time).

- count, population:

  Columns of counts and populations.

- age:

  Column of age groups. For clustering their order matters: the levels
  of a factor, numeric order, or the leading number of each label (for
  example `"0-4"`, `"5-9"`, `"85+"`).

- standard:

  Either a column of `data` holding the standard population of each row,
  or a numeric vector named by age group (for example from
  [std_populations](https://skhademioore.github.io/joinpoint/reference/std_populations.md)).

- by:

  Optional character vector of by-variables.

- delay:

  Optional column of delay factors (1 to 20).

- rate_per:

  Rate multiplier.

## Value

A data frame with the by-variables, `x`, `rate`, `se`, `count` and
`population`, ready for
[`joinpoint()`](https://skhademioore.github.io/joinpoint/reference/joinpoint.md)
with `y = rate, se = se`.

## Details

For age groups `a` with counts `C_a`, populations `P_a`, delay factors
`D_a` and standard population weights `w_a = S_a / sum(S)`, the rate is
`M * sum(w_a * D_a * C_a / P_a)` with variance
`M^2 * sum(w_a^2 * D_a^2 * C_a / P_a^2)`, assuming Poisson counts.

## Examples

``` r
us <- subset(std_populations, standard == "segi1960")
std <- stats::setNames(us$population, us$age_group)
d <- expand.grid(year = 2000:2004, age = names(std))
d$count <- 10 + seq_len(nrow(d)) %% 7
d$pop <- 1e5
age_adjust(d, year, count, pop, age, std)
#>   year   rate        se count population
#> 1 2000 13.000 0.9857738   235    1800000
#> 2 2001 12.810 0.9750897   232    1800000
#> 3 2002 12.970 0.9847589   236    1800000
#> 4 2003 12.885 0.9811345   233    1800000
#> 5 2004 13.255 1.0030329   237    1800000
```
