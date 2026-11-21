## Submission

This is a new submission of joinpoint 0.1.0.

## Test environments

* GitHub Actions: Ubuntu (R-devel, R-release, R-oldrel-1), macOS (R-release)
  and Windows (R-release).
* Local: macOS 26.5, aarch64, R 4.5.1, with the CRAN incoming checks.

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission.

`R CMD check --as-cran` gave 0 errors, 0 warnings and 0 notes on every
GitHub Actions platform. The local run with the CRAN incoming checks gave
the note above. Locally, Apple clang 21 also warned about R's own header
`R_ext/Boolean.h` (unknown warning group `-Wfixed-enum-extension`), which
is unrelated to the package.

Possibly misspelled words in DESCRIPTION, if reported, are terms of the
method (joinpoint, joinpoints, AAPC, autocorrelated) and the names of the
cited authors (Lerman, Czajkowski, et al.).

## Notes for the reviewers

* Random numbers: the permutation test and the empirical quantile intervals
  are reproducible by default through the documented arguments
  `jp_control(perm_seed = 7160, eq_seed = 10000)`. The functions save the
  global random number state before setting a seed and restore it on exit,
  so the user's random number stream is not changed (except for the normal
  deviate cached by the Box-Muller generator, which R does not expose);
  `perm_seed = NULL` or `eq_seed = NULL` uses the current stream instead.
  Interactive graphics do not use the random number generator.
* Run time: examples take at most about one second each. Slow permutation
  tests are skipped on CRAN; the full test suite and the vignettes take
  about two minutes of CPU time in total.
* Files are written only by `jp_export()`, `jp_save()` and `jp_report()`,
  to paths the user supplies; the examples, tests and vignettes write to
  `tempdir()` and remove what they write. Documents are rendered in a
  private temporary folder and then copied to the requested path, so no
  other file is replaced or removed.
* Rendered reports (`jp_report()`) and document exports use Pandoc through
  the suggested rmarkdown package, only when it is available; reports are
  rendered by a separate `Rscript` process (without user or site profiles),
  so that a report can be made from within a document being knitted. The
  `\donttest` report example takes about five seconds.
* Interactive graphics use the suggested packages ggiraph and gdtools when
  installed; they name system fonts and do not write to the gdtools font
  cache.
* References for the methods are given in the Description field with DOIs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
