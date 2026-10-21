# Control settings for joinpoint fitting

Settings of the grid search, model selection, confidence intervals, AAPC
ranges and random number generation.

## Usage

``` r
jp_control(
  min_obs_end = 2,
  min_obs_between = 2,
  n_grid = 0,
  alpha = 0.05,
  n_perm = 4499,
  n_exceed = NULL,
  perm_seed = 7160,
  n_resample = 5001,
  eq_seed = 10000,
  eq_final_only = TRUE,
  seed_mode = c("constant", "vary"),
  alpha_apc = 0.05,
  alpha_aapc = 0.05,
  alpha_joinpoint = 0.05,
  alpha_jump = 0.05,
  aapc_full = TRUE,
  aapc_last = NULL,
  aapc_ranges = NULL,
  madwd = NULL,
  dds_cutoffs = c(10, 200),
  infeasible = c("reduce", "recommended", "skip"),
  tau_interval = c("lerman", "literal"),
  time_budget = 60
)
```

## Arguments

- min_obs_end:

  Minimum number of observations from a joinpoint to either end of the
  data, excluding the joinpoint itself when it falls on an observation
  (1 to 9, default 2).

- min_obs_between:

  Minimum number of observations between two joinpoints, excluding
  joinpoints that fall on observations (0 to 9, default 2).

- n_grid:

  Number of grid points placed between adjacent observed x values in the
  grid search (0 to 9, default 0).

- alpha:

  Overall significance level of the sequential permutation test.

- n_perm:

  Number of randomly permuted data sets for the permutation test
  (default 4499).

- n_exceed:

  Optional early stopping of the permutation and likelihood ratio tests
  (Besag and Clifford 1991): simulation stops once `n_exceed` simulated
  statistics reach the observed one, at simulation `L`, and the p value
  is `n_exceed / L`. `NULL` runs all `n_perm` simulations. Values of 10
  to 20 save most of the time of tests that do not reject.

- perm_seed:

  Random number seed for the permutation test (`NULL` for the current
  random number stream). The global random number state is saved before
  the seed is set and restored afterwards.

- n_resample:

  Number of resampled data sets for the empirical quantile confidence
  intervals (default 5001).

- eq_seed:

  Random number seed for the empirical quantile method (`NULL` for the
  current random number stream). The global random number state is saved
  before the seed is set and restored afterwards.

- eq_final_only:

  Use the empirical quantile method only for the final selected model;
  other models get parametric intervals.

- seed_mode:

  `"constant"` uses the same seed for every cohort; `"vary"` uses
  `seed * g + (g - 1)` (modulo `2^31 - 1`) for cohort `g`.

- alpha_apc, alpha_aapc, alpha_joinpoint, alpha_jump:

  Significance levels for the APC, AAPC, joinpoint location and jump or
  comparability ratio intervals.

- aapc_full:

  Compute the AAPC over the full range of the data.

- aapc_last:

  Integer vector (up to two values) of "Last X observations" AAPC
  ranges.

- aapc_ranges:

  List of up to seven numeric `c(start, end)` ranges anywhere within the
  data range.

- madwd:

  Minimum APC difference worth detecting, in percentage points (`NULL`
  to disable).

- dds_cutoffs:

  Cutoffs `c < d` of the data dependent selection.

- infeasible:

  What to do when a cohort has too few observations for the requested
  `kmax` under the placement constraints: `"reduce"` lowers `kmax` to
  the largest feasible value, `"recommended"` uses the recommended
  maximum for the cohort, and `"skip"` leaves the cohort out. A warning
  is given in every case.

- tau_interval:

  Parametric interval of joinpoint locations: `"lerman"` (Lerman 1980)
  or the more conservative `"literal"` variant `Q (1 + k / (n - p)) F`,
  which reproduces some published intervals.

- time_budget:

  Approximate computing time, in seconds, allowed for the exhaustive
  search when `kmax` is left to its recommended value. The recommended
  maximum is lowered, with a message, when the search (with its
  permutation tests or resamples) would take longer; the budget is an
  estimate, not a deadline, and `Inf` removes it. An explicit `kmax` is
  always honored.

## Value

A list of class `jp_control`.

## References

Besag J, Clifford P (1991). Sequential Monte Carlo p-values. Biometrika
78:301-304.

## Examples

``` r
jp_control(n_grid = 1, aapc_last = c(5, 10))
#> $min_obs_end
#> [1] 2
#> 
#> $min_obs_between
#> [1] 2
#> 
#> $n_grid
#> [1] 1
#> 
#> $alpha
#> [1] 0.05
#> 
#> $n_perm
#> [1] 4499
#> 
#> $n_exceed
#> NULL
#> 
#> $perm_seed
#> [1] 7160
#> 
#> $n_resample
#> [1] 5001
#> 
#> $eq_seed
#> [1] 10000
#> 
#> $eq_final_only
#> [1] TRUE
#> 
#> $seed_mode
#> [1] "constant"
#> 
#> $alpha_apc
#> [1] 0.05
#> 
#> $alpha_aapc
#> [1] 0.05
#> 
#> $alpha_joinpoint
#> [1] 0.05
#> 
#> $alpha_jump
#> [1] 0.05
#> 
#> $aapc_full
#> [1] TRUE
#> 
#> $aapc_last
#> [1]  5 10
#> 
#> $aapc_ranges
#> NULL
#> 
#> $madwd
#> NULL
#> 
#> $dds_cutoffs
#> [1]  10 200
#> 
#> $infeasible
#> [1] "reduce"
#> 
#> $tau_interval
#> [1] "lerman"
#> 
#> $time_budget
#> [1] 60
#> 
#> attr(,"class")
#> [1] "jp_control"
```
