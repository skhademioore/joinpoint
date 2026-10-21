## Submission

This is a new submission of joinpoint 0.1.0.

## R CMD check results

`R CMD check --as-cran` on the built tarball (macOS, aarch64, R 4.6.0) gives
0 errors, 0 warnings and 1 note:

* New submission.
* Possibly misspelled words in DESCRIPTION: these are the method's name
  (joinpoint, joinpoints, AAPC, autocorrelated) and the authors of the cited
  papers (Lerman, Czajkowski, et al.).

The GitHub URLs in DESCRIPTION, `inst/CITATION` and the package help page
resolve once the repository is public.

## Notes for the reviewers

* Random numbers: the permutation test and the empirical quantile intervals
  are reproducible by default through the documented arguments
  `jp_control(perm_seed = 7160, eq_seed = 10000)`. The functions save the
  global random number state before setting a seed and restore it on exit,
  so the user's random number stream is not changed; `perm_seed = NULL` or
  `eq_seed = NULL` uses the current stream instead.
* Run time: examples take at most about one second each. Slow permutation
  tests are skipped on CRAN; the full test suite and the vignettes take
  about two minutes of CPU time in total.
* Files are written only by `jp_export()`, to a directory the user must
  supply; the example and the tests write to `tempdir()` and remove what they
  write.
* References for the methods are given in the Description field with DOIs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
