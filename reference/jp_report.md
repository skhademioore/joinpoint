# Joinpoint report as a document, R Markdown or Quarto file

Writes a report of a joinpoint fit: the settings of the analysis, the
trends of the selected models (APCs and AAPCs with confidence
intervals), the joinpoints, the model selection statistics and the
figures. The format follows the extension of `file`:

## Usage

``` r
jp_report(
  object,
  file,
  title = "Joinpoint regression analysis",
  author = NULL,
  model = "final",
  plots = c("trend", "apc", "aapc", "selection"),
  format = c("html", "docx", "pdf"),
  interactive = NULL,
  quiet = TRUE
)
```

## Arguments

- object:

  A `joinpoint` object.

- file:

  Output file; the extension sets the format.

- title, author:

  Title and author of the report.

- model:

  `"final"` or a number of joinpoints.

- plots:

  Figures to include: any of `"trend"`, `"apc"`, `"aapc"` and
  `"selection"`, or none
  ([`character()`](https://rdrr.io/r/base/character.html)).

- format:

  Output format written into the header of an `.Rmd` or `.qmd` source:
  `"html"`, `"docx"` or `"pdf"`.

- interactive:

  Interactive figures (HTML output only); by default when the output is
  HTML and the ggiraph package is installed.

- quiet:

  Render without printing the progress of knitr and Pandoc.

## Value

The path of the written file, invisibly.

## Details

- `.Rmd` or `.qmd`: an R Markdown or Quarto source document, with the
  fit saved beside it (`<name>_fit.rds`) and read by its first chunk,
  ready to edit and render (in RStudio, with
  [`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html),
  or with `quarto render`).

- `.html`, `.docx`, `.odt`, `.pdf`, `.md` or `.tex`: the rendered report
  (with the rmarkdown package and Pandoc; PDF also needs LaTeX). The
  figures of HTML reports are interactive when the ggiraph package is
  installed; the LaTeX and Markdown reports keep their figures in a
  folder beside them.

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
dir <- file.path(tempdir(), "joinpoint-report")
dir.create(dir)
jp_report(fit, file.path(dir, "tuberculosis.qmd"))
list.files(dir)
#> [1] "tuberculosis.qmd"     "tuberculosis_fit.rds"
# \donttest{
if (requireNamespace("rmarkdown", quietly = TRUE) &&
    rmarkdown::pandoc_available()) {
  jp_report(fit, file.path(dir, "tuberculosis.docx"))
}
# }
unlink(dir, recursive = TRUE)
```
