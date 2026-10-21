# Standard populations for age adjustment

Standard population distributions by age group for direct age adjustment
with
[`age_adjust()`](https://skhademioore.github.io/joinpoint/reference/age_adjust.md)
and
[`jp_cluster()`](https://skhademioore.github.io/joinpoint/reference/jp_cluster.md).

## Usage

``` r
std_populations
```

## Format

A data frame with columns `standard`, `age_group` and `population`.
Standards: `"us2000"` (2000 United States standard population, Census
P25-1130, 19 age groups), `"segi1960"` (Segi world standard, 18 groups,
per 100,000) and `"esp2013"` (2013 European standard population, 19
groups, per 100,000).

## Source

United States Census Bureau, Current Population Reports P25-1130 (2000
standard); Segi M (1960), Cancer mortality for selected sites in 24
countries; Eurostat (2013), Revision of the European Standard
Population.

## Examples

``` r
head(std_populations)
#>   standard age_group population
#> 1   us2000        00    3794901
#> 2   us2000     01-04   15191619
#> 3   us2000     05-09   19919840
#> 4   us2000     10-14   20056779
#> 5   us2000     15-19   19819518
#> 6   us2000     20-24   18257225
tapply(std_populations$population, std_populations$standard, sum)
#>   esp2013  segi1960    us2000 
#>    100000    100000 274633642 
```
