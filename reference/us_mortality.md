# Age adjusted death rates for six causes, United States 1900 to 1998

Age adjusted death rates per 100,000 for accidents, cancer, heart
disease, influenza and pneumonia, stroke and tuberculosis in the Death
Registration States (1900 to 1932) and the United States (1933 to 1998),
with the resident population.

## Usage

``` r
us_mortality
```

## Format

A data frame with 594 rows and 4 columns:

- year:

  Calendar year.

- cause:

  Cause of death (factor).

- rate:

  Age adjusted death rate per 100,000.

- population:

  United States resident population.

## Source

CDC National Center for Health Statistics, National Vital Statistics
System, Age Adjusted Death Rates for Selected Causes, Death Registration
States 1900 to 1932 and United States 1933 to 1998 (HIST293), as
compiled by van Doren T (2024), Estimating time points of significant
change in cause-specific mortality: joinpoint regression in R, The
Denominator, [doi:10.6069/JK1Y-7Z16](https://doi.org/10.6069/JK1Y-7Z16)
.

## Examples

``` r
head(us_mortality)
#>   year     cause  rate population
#> 1 1900 Accidents  90.3   76090000
#> 2 1901 Accidents 109.3   77580000
#> 3 1902 Accidents  93.6   79160000
#> 4 1903 Accidents 106.9   80630000
#> 5 1904 Accidents 112.8   82170000
#> 6 1905 Accidents 108.5   83820000
table(us_mortality$cause)
#> 
#>               Accidents                  Cancer           Heart Disease 
#>                      99                      99                      99 
#> Influenza and Pneumonia                  Stroke            Tuberculosis 
#>                      99                      99                      99 
```
