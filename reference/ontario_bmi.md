# Prepregnancy body mass index in Ontario, 2012 to 2021

Prepregnancy body mass index (BMI) of 1,170,997 singleton births in
Ontario, Canada, by fiscal year (April to March) from 2012-2013 to
2021-2022: the number of births in each World Health Organization BMI
category, and the yearly mean and standard deviation of BMI.

## Usage

``` r
ontario_bmi
```

## Format

A data frame with 60 rows and 6 columns:

- year:

  Start of the fiscal year (2012 is April 2012 to March 2013).

- category:

  BMI category (ordered from underweight to obesity class III).

- count:

  Births in the category.

- births:

  All births in the fiscal year.

- mean_bmi,sd_bmi:

  Yearly mean and standard deviation of BMI in kg/m^2, rounded to one
  decimal and repeated over categories.

## Source

Khademioore S, Darling EK, D'Souza R, Muraca GM, Guo Y, Anderson LN
(2026). Prepregnancy body mass index in Ontario from 2012 to 2022: a
population-based descriptive study. Canadian Journal of Public Health.
[doi:10.17269/s41997-026-01214-w](https://doi.org/10.17269/s41997-026-01214-w)
. Table 3, licensed under CC BY 4.0.

## Examples

``` r
head(ontario_bmi)
#>   year          category count births mean_bmi sd_bmi
#> 1 2012       Underweight  9262 107457     25.3    6.2
#> 2 2012     Normal weight 53904 107457     25.3    6.2
#> 3 2012        Overweight 21955 107457     25.3    6.2
#> 4 2012   Obesity class I 14541 107457     25.3    6.2
#> 5 2012  Obesity class II  3957 107457     25.3    6.2
#> 6 2012 Obesity class III  3838 107457     25.3    6.2
```
