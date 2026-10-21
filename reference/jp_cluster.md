# Multi-group clustering of age specific trends

Clusters contiguous age groups whose trends are as parallel as possible
and fits a joinpoint model to the age adjusted rates of each cluster:
pre-analysis count compression, parallel joinpoint fits for every block
of contiguous age groups, the partition into `m` clusters with minimum
total residual sum of squares, selection of `m` by BIC, and optional
clustering MADWD.

## Usage

``` r
jp_cluster(
  data,
  x,
  count,
  population,
  age,
  standard,
  by = NULL,
  delay = NULL,
  min_count = 0,
  compression = c("forward", "backward"),
  method = c("wbic", "bic", "bic3", "permutation"),
  max_clusters = NULL,
  madwd = NULL,
  kmin = 0,
  kmax = NULL,
  final_method = "wbic",
  rate_per = 1e+05,
  control = jp_control(),
  ...
)
```

## Arguments

- data:

  A data frame with one row per age group and time point (and cohort).

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

- min_count:

  Age groups are combined until counts exceed this value at every time
  point.

- compression:

  `"forward"` or `"backward"` count compression.

- method:

  Method for the number of joinpoints within clustering: `"wbic"`,
  `"bic"`, `"bic3"` or `"permutation"`.

- max_clusters:

  Maximum number of clusters (default: number of compressed age groups).

- madwd:

  Clustering minimum APC difference worth detecting, in percentage
  points (`NULL` to disable).

- kmin, kmax:

  Minimum and maximum number of joinpoints.

- final_method:

  Model selection method for the final fits of the cluster age adjusted
  rates (see
  [`joinpoint()`](https://skhademioore.github.io/joinpoint/reference/joinpoint.md)).

- rate_per:

  Rate multiplier.

- control:

  Settings from
  [`jp_control()`](https://skhademioore.github.io/joinpoint/reference/jp_control.md),
  used for the clustering fits and the final fits.

- ...:

  Further arguments to
  [`joinpoint()`](https://skhademioore.github.io/joinpoint/reference/joinpoint.md)
  for the final cluster fits (for example `ci_method`).

## Value

An object of class `jp_cluster` with the compressed age groups, the
cluster selection table, the final clusters, the MADWD steps, the
cluster rates and `fit`, the
[`joinpoint()`](https://skhademioore.github.io/joinpoint/reference/joinpoint.md)
fit of the cluster age adjusted rates.

## Examples

``` r
set.seed(1)
ages <- c("40-49", "50-59", "60-69", "70-79")
d <- expand.grid(year = 2000:2019, age = ages, stringsAsFactors = FALSE)
d$pop <- 1e5
slope <- ifelse(d$age %in% ages[1:2], 0.03, -0.02)
d$count <- rpois(nrow(d), 200 * exp(slope * (d$year - 2000)))
cl <- jp_cluster(d, year, count, pop, age, standard = setNames(rep(1, 4), ages),
                 kmax = 1, ci_method = "parametric")
cl$clusters
#>          cluster age_group
#> 1 40-49 to 50-59     40-49
#> 2 40-49 to 50-59     50-59
#> 3 60-69 to 70-79     60-69
#> 4 60-69 to 70-79     70-79
```
