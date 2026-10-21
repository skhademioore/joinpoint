# Contributing to joinpoint

Thank you for considering a contribution to `joinpoint`. The package implements
statistical methods for trend analysis, so changes need both software-quality
and statistical-methodology review.

## Before You Start

- Use GitHub issues for bug reports, feature requests, documentation gaps, and
  methodological questions.
- For larger changes, open an issue first so the design can be discussed before
  implementation.
- Keep pull requests focused. A statistical-method change, a documentation
  update, and a refactor should usually be separate pull requests.

## Development Setup

Install the package dependencies from the repository root:

```r
install.packages(c("remotes", "testthat", "devtools", "roxygen2"))
remotes::install_deps(dependencies = TRUE)
```

The grid search is compiled C++ (through Rcpp), so you need a working C++
toolchain:

- macOS: Xcode Command Line Tools (`xcode-select --install`).
- Windows: Rtools.
- Linux: `g++` and `make`.

## Branches and Pull Requests

- One branch per logical change, with a name that describes it
  (`aapc-ratio-comparison`, `continuous-interval-limits`).
- Write commit subjects in the imperative ("Add the jump model"). Conventional
  Commit prefixes are not used here.
- Never rewrite `main`. `git commit --amend` and `git push --force-with-lease`
  are acceptable only on an unmerged branch that nobody else has pulled.
- Keep the implementation, its tests, the roxygen source, the regenerated
  documentation, and the `NEWS.md` entry together in the same pull request.

## Code Style

- Follow the style of the surrounding R files.
- Prefer clear validation and explicit errors over silent coercion.
- Do not add package dependencies unless they are necessary and justified.
- Every function and line should have a reason to exist; prefer a small,
  correct implementation over a general one.
- Write documentation in American English, and cite academic papers for
  methods.
- Update roxygen comments and regenerate the `.Rd` files together.

## Testing

Run the tests while developing:

```r
devtools::test()
```

Before opening a pull request, run the full check:

```sh
R CMD build .
R CMD check --as-cran joinpoint_*.tar.gz
```

A change to an estimator, a test, or an interval needs a test that checks the
statistical output against an independent computation (brute force, `lm()`,
`glm()`, or a published analysis), not only that the code runs.

## Documentation

Update user-facing documentation when behavior changes:

- roxygen examples and parameter descriptions,
- README examples,
- vignettes, and their rendered HTML in `vignettes/` (see below),
- `NEWS.md`.

Every number quoted in a vignette must match the output it accompanies;
prefer inline R code for numbers in prose. To refresh the rendered vignettes:

```r
for (f in list.files("vignettes", "\\.Rmd$", full.names = TRUE)) {
  rmarkdown::render(f, quiet = TRUE)
}
```

## Pull Request Checklist

Before opening a pull request, confirm:

- Tests relevant to the change pass, and `R CMD check` is clean.
- Documentation is updated.
- `NEWS.md` has an entry for user-facing changes.
- No build artifacts, `.Rcheck` directories, compiled objects, or tarballs are
  committed.
- The description explains both the software change and, when relevant, its
  statistical impact.

## Code of Conduct

All contributors are expected to follow the project
[Code of Conduct](CODE_OF_CONDUCT.md).
