# Export joinpoint results

Writes the result tables of a joinpoint fit to delimited text files, one
file per table: data and modeled values, APCs, AAPCs, joinpoint
locations, model estimates, model selection (criteria, permutation
tests, data dependent selection, MADWD and the final selected models
with their search limits), coefficient covariances, jump estimates,
pairwise and AAPC comparisons, clustering tables, skipped cohorts, the
centered covariate design, a supplied covariance matrix and every
setting of the run. [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)
keeps the complete fitted object.

## Usage

``` r
jp_export(
  object,
  dir,
  prefix = "joinpoint",
  model = c("all", "final"),
  comparison = NULL,
  sep = ","
)
```

## Arguments

- object:

  A `joinpoint` or `jp_cluster` object.

- dir:

  Output directory (created if needed).

- prefix:

  File name prefix.

- model:

  `"all"` (every fitted model) or `"final"`.

- comparison:

  Optional result of
  [`jp_compare()`](https://skhademioore.github.io/joinpoint/reference/jp_compare.md)
  to export with the fit.

- sep:

  Field separator.

## Value

The paths of the written files, invisibly.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
dir <- file.path(tempdir(), "joinpoint-example")
files <- jp_export(fit, dir)
basename(files)
#> [1] "joinpoint.Data.txt"                  
#> [2] "joinpoint.APCs.txt"                  
#> [3] "joinpoint.AAPCs.txt"                 
#> [4] "joinpoint.Joinpoints.txt"            
#> [5] "joinpoint.Model.Estimates.txt"       
#> [6] "joinpoint.Coefficient.Covariance.txt"
#> [7] "joinpoint.Model.Selection.txt"       
#> [8] "joinpoint.Final.Selected.Model.txt"  
#> [9] "joinpoint.Run.Summary.txt"           
unlink(dir, recursive = TRUE)
```
