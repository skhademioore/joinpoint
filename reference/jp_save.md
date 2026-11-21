# Save a table, figure or report

Writes a result table, a figure or an interactive graphic to a file
whose format follows its extension, for a manuscript, a report or a
presentation.

## Usage

``` r
jp_save(
  x,
  file,
  title = NULL,
  note = NULL,
  digits = 3,
  width = 7.5,
  height = 5,
  dpi = 300,
  label = NULL,
  ...
)
```

## Arguments

- x:

  A data frame, a ggplot, an htmlwidget or a `joinpoint` fit.

- file:

  Output file; the extension sets the format.

- title:

  Optional caption of the table or figure.

- note:

  Optional note printed under the table or figure.

- digits:

  Decimals of the numeric columns of a table (character columns, such as
  those of
  [`jp_table()`](https://skhademioore.github.io/joinpoint/reference/jp_table.md),
  are written as they are). In `.csv`, `.tsv` and `.txt` files numbers
  stay numeric (p values rounded to `digits` significant digits); other
  formats show p values with three decimals, and below 0.001 as "\<
  0.001". Line breaks inside cells of Markdown and document tables
  become spaces.

- width, height, dpi:

  Size of a figure in inches and its resolution.

- label:

  Optional LaTeX label (`.tex` output).

- ...:

  Further arguments to
  [`jp_report()`](https://skhademioore.github.io/joinpoint/reference/jp_report.md)
  when `x` is a fit.

## Value

The path of the written file, invisibly.

## Details

- Tables (data frames such as the results of
  [`jp_table()`](https://skhademioore.github.io/joinpoint/reference/jp_table.md),
  [`apc()`](https://skhademioore.github.io/joinpoint/reference/apc.md)
  or
  [`model_selection()`](https://skhademioore.github.io/joinpoint/reference/model_selection.md)`$criteria`):
  `.csv`, `.tsv`, `.txt`, `.md` (Markdown pipe table), `.tex` (LaTeX
  `table` with booktabs rules; load `\usepackage{booktabs}`), `.html`,
  `.docx` (Word table), `.odt`, `.rtf`, `.pptx` and `.pdf`.

- Figures (ggplot objects such as
  [autoplot()](https://skhademioore.github.io/joinpoint/reference/autoplot.joinpoint.md)
  or
  [`plot_apc()`](https://skhademioore.github.io/joinpoint/reference/plot_apc.md)):
  the image formats of
  [`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
  (`.png`, `.pdf`, `.svg`, `.jpeg`, `.tiff`, `.eps`, ...), `.tex` (a
  LaTeX `figure` that includes a PDF of the figure written beside it),
  and `.html`, `.docx`, `.odt`, `.rtf`, `.pptx` and `.md` documents
  holding the figure.

- Interactive graphics
  ([`jp_interactive()`](https://skhademioore.github.io/joinpoint/reference/jp_interactive.md)):
  `.html`, a standalone web page (without Pandoc, a page with its
  scripts in the folder `<name>_files` beside it).

- A `joinpoint` fit: a complete report, see
  [`jp_report()`](https://skhademioore.github.io/joinpoint/reference/jp_report.md).

Word, OpenDocument, PowerPoint, RTF, HTML and PDF documents are made
with Pandoc, which comes with RStudio and Quarto (or see
<https://pandoc.org>), through the rmarkdown package; PDF documents also
need LaTeX (for example
[`tinytex::install_tinytex()`](https://rdrr.io/pkg/tinytex/man/install_tinytex.html)).

## Examples

``` r
tb <- subset(us_mortality, cause == "Tuberculosis")
fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
dir <- file.path(tempdir(), "joinpoint-save")
dir.create(dir)
jp_save(jp_table(fit), file.path(dir, "trends.tex"),
        title = "Trends in tuberculosis mortality")
jp_save(jp_table(fit), file.path(dir, "trends.md"))
jp_save(autoplot(fit), file.path(dir, "trend.png"), width = 6, height = 4)
if (requireNamespace("rmarkdown", quietly = TRUE) &&
    rmarkdown::pandoc_available()) {
  jp_save(jp_table(fit), file.path(dir, "trends.docx"))
}
list.files(dir)
#> [1] "trend.png"   "trends.docx" "trends.md"   "trends.tex" 
unlink(dir, recursive = TRUE)
```
