# Sharing results: tables, figures, reports and interactive graphics

An analysis ends in a manuscript, a report, a slide or a web page. This
vignette shows how the tables and figures of a joinpoint fit get there:
saved to files in the format each needs, written into a complete report,
printed in R Markdown and Quarto documents, and drawn as interactive
graphics.

The examples use the tuberculosis and stroke mortality of
`us_mortality`:

``` r

library(joinpoint)
fit <- joinpoint(subset(us_mortality, cause %in% c("Stroke", "Tuberculosis")),
                 year, rate, by = "cause", kmax = 3, ci_method = "parametric")
dir <- file.path(tempdir(), "joinpoint-reports")
dir.create(dir, showWarnings = FALSE)
```

## Tables

[`jp_save()`](https://skhademioore.github.io/joinpoint/reference/jp_save.md)
writes a table in the format given by the extension of the file: `.csv`,
`.tsv` and `.txt` for spreadsheets and other software, `.md` for
Markdown, `.tex` for LaTeX, and `.docx`, `.odt`, `.rtf`, `.pptx`,
`.html` and `.pdf` for documents. A title becomes the caption, and a
note is printed under the table:

``` r

tab <- jp_table(fit)
jp_save(tab, file.path(dir, "trends.tex"), title = "Trends in mortality",
        note = "* The APC or AAPC differs from 0 at alpha = 0.05.",
        label = "tab:trends")
```

The LaTeX table uses the rules of the booktabs package
(`\usepackage{booktabs}` in the preamble) and escapes the special
characters of LaTeX:

    \begin{table}[htbp]
    \centering
    \caption{Trends in mortality}
    \label{tab:trends}
    \begin{tabular}{llrrr}
    \toprule
    cause & Trend & Period & APC / AAPC (95\% CI) & p \\
    \midrule
    Stroke & Segment 1 & 1900-1945 & -0.94* (-1.06, -0.83) & $<$ 0.001 \\
    Stroke & Segment 2 & 1945-1968 & -0.05 (-0.40, 0.30) & 0.757 \\
    Stroke & Segment 3 & 1968-1987 & -4.68* (-5.13, -4.23) & $<$ 0.001 \\
    Stroke & Segment 4 & 1987-1998 & -1.29* (-2.27, -0.30) & 0.011 \\
    Stroke & AAPC & 1900-1998 & -1.51* (-1.68, -1.34) & $<$ 0.001 \\
    Tuberculosis & Segment 1 & 1900-1947 & -3.72* (-3.91, -3.53) & $<$ 0.001 \\
    Tuberculosis & Segment 2 & 1947-1954 & -14.66* (-18.35, -10.80) & $<$ 0.001 \\
    Tuberculosis & Segment 3 & 1954-1981 & -8.97* (-9.41, -8.53) & $<$ 0.001 \\
    Tuberculosis & Segment 4 & 1981-1998 & -4.52* (-5.39, -3.64) & $<$ 0.001 \\
    Tuberculosis & AAPC & 1900-1998 & -6.14* (-6.50, -5.78) & $<$ 0.001 \\
    \bottomrule
    \end{tabular}
    \par\smallskip{\footnotesize * The APC or AAPC differs from 0 at alpha = 0.05.}
    \end{table}

Any data frame can be saved this way, for example the APCs as numbers or
the model selection statistics, with the decimals set by `digits`:

``` r

jp_save(apc(fit), file.path(dir, "apc.csv"))
jp_save(model_selection(fit)$criteria, file.path(dir, "selection.md"),
        digits = 2)
```

Word, OpenDocument, PowerPoint, RTF, HTML and PDF documents are made
with Pandoc, which comes with RStudio and Quarto, through the rmarkdown
package; PDF documents also need LaTeX (for example
[`tinytex::install_tinytex()`](https://rdrr.io/pkg/tinytex/man/install_tinytex.html)).
A Word table for a manuscript:

``` r

jp_save(tab, file.path(dir, "trends.docx"), title = "Trends in mortality")
```

## Figures

Figures are saved the same way: in the image formats of
[`ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html)
(`.png`, `.pdf`, `.svg`, `.tiff`, `.eps` and others), as a LaTeX figure,
or inside a Word, OpenDocument, PowerPoint or HTML document. `width`,
`height` (inches) and `dpi` set the size and resolution:

``` r

g <- autoplot(fit, color_by = "trend", log_y = TRUE)
jp_save(g, file.path(dir, "trends.png"), width = 7, height = 4, dpi = 300)
jp_save(g, file.path(dir, "trends.tex"), title = "Mortality trends",
        label = "fig:trends")
```

The `.tex` file is a `figure` environment that includes a PDF of the
figure written beside it (`trends.pdf`), ready for `\input{}` in a
manuscript.

## Complete reports

[`jp_report()`](https://skhademioore.github.io/joinpoint/reference/jp_report.md)
writes the whole analysis: the settings in words, the table of APCs and
AAPCs, the joinpoints, the model selection statistics and the figures.
An `.Rmd` or `.qmd` file gives the R Markdown or Quarto source, with the
fit saved beside it, ready to edit and render:

``` r

jp_report(fit, file.path(dir, "mortality.qmd"),
          title = "Trends in mortality", author = "A. Analyst")
```

    ---
    title: "Trends in mortality"
    author: "A. Analyst"
    date: "2026-11-21"
    format:
      html:
        embed-resources: true
    execute:
      echo: false
      warning: false
      message: false
    ---

    ```{r setup, include = FALSE}
    library(joinpoint)
    fit <- readRDS("mortality_fit.rds")
    ```

    ## Analysis

    The trend of Rate over Year was modeled by a log-linear joinpoint model with constant variance for each of 2 cohorts. The number of joinpoints (0 to at most 3) was selected by the weighted BIC; the selected models have 3 joinpoints. The tables give parametric confidence intervals at the 95% level, and the 95% level for the joinpoints. The intervals are conditional on the number of joinpoints.

    ## Annual Percent Changes

    ```{r}
    knitr::kable(jp_table(fit), caption = "APCs of the segments and AAPCs, with 95% confidence intervals.")
    ```

    ## Joinpoints

    ```{r}
    jt <- joinpoints(fit)

Other extensions give the rendered report: `.html`, `.docx`, `.odt`,
`.pdf`, `.md` or `.tex`. The figures of HTML reports are interactive
when the ggiraph package is installed.

``` r

jp_report(fit, file.path(dir, "mortality.docx"), title = "Trends in mortality")
```

[`jp_save()`](https://skhademioore.github.io/joinpoint/reference/jp_save.md)
of a fit writes the same report.

## Fits in R Markdown and Quarto documents

A fit printed in a chunk of an R Markdown or Quarto document shows the
same summary as in the console. With the option `joinpoint.knit_table`,
set once in the first chunk of the document, it appears instead as the
table of
[`jp_table()`](https://skhademioore.github.io/joinpoint/reference/jp_table.md),
formatted for the output (HTML, LaTeX, Word or Markdown); the chunk
option `tab.cap` sets its caption:

``` r

options(joinpoint.knit_table = TRUE)
fit
```

| cause        | Trend     | Period    | APC / AAPC (95% CI)       | p        |
|:-------------|:----------|:----------|:--------------------------|:---------|
| Stroke       | Segment 1 | 1900-1945 | -0.94\* (-1.06, -0.83)    | \< 0.001 |
| Stroke       | Segment 2 | 1945-1968 | -0.05 (-0.40, 0.30)       | 0.757    |
| Stroke       | Segment 3 | 1968-1987 | -4.68\* (-5.13, -4.23)    | \< 0.001 |
| Stroke       | Segment 4 | 1987-1998 | -1.29\* (-2.27, -0.30)    | 0.011    |
| Stroke       | AAPC      | 1900-1998 | -1.51\* (-1.68, -1.34)    | \< 0.001 |
| Tuberculosis | Segment 1 | 1900-1947 | -3.72\* (-3.91, -3.53)    | \< 0.001 |
| Tuberculosis | Segment 2 | 1947-1954 | -14.66\* (-18.35, -10.80) | \< 0.001 |
| Tuberculosis | Segment 3 | 1954-1981 | -8.97\* (-9.41, -8.53)    | \< 0.001 |
| Tuberculosis | Segment 4 | 1981-1998 | -4.52\* (-5.39, -3.64)    | \< 0.001 |
| Tuberculosis | AAPC      | 1900-1998 | -6.14\* (-6.50, -5.78)    | \< 0.001 |

Trends in stroke and tuberculosis mortality. {.table}

Figures are drawn by
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
and the other plot functions as in any chunk.

## Interactive graphics

The trend figure
([`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)),
the APC and AAPC figures
([`plot_apc()`](https://skhademioore.github.io/joinpoint/reference/plot_apc.md),
[`plot_aapc()`](https://skhademioore.github.io/joinpoint/reference/plot_apc.md))
and the model selection figure
([`plot_selection()`](https://skhademioore.github.io/joinpoint/reference/plot_selection.md))
can be drawn as interactive graphics, with `interactive = TRUE` in the
plot function and
[`jp_interactive()`](https://skhademioore.github.io/joinpoint/reference/jp_interactive.md),
or
[`jp_interactive()`](https://skhademioore.github.io/joinpoint/reference/jp_interactive.md)
on the fit directly, with `plot` choosing the figure. The graphic looks
like the static figure. Hovering a segment shows its period, APC,
confidence interval and p value and highlights the segment with its band
and label; hovering an observation shows its value, its confidence
interval and the fitted trend; hovering a joinpoint shows its location
and confidence interval. The toolbar zooms and saves the graphic as an
image.

``` r

jp_interactive(fit, color_by = "trend", log_y = TRUE, height = 4.2)
```

The other figures are chosen with `plot`, here the APC of each segment:

``` r

jp_interactive(fit, plot = "apc", height = 3.5)
```

A figure changed with the grammar of ggplot2 stays interactive:

``` r

g <- autoplot(fit, interactive = TRUE) +
  ggplot2::labs(title = "Mortality per 100,000", x = NULL)
jp_interactive(g, height = 4.2)
```

The graphics are htmlwidgets: they appear in the RStudio viewer, in HTML
documents made with R Markdown or Quarto, on pkgdown sites and in Shiny
applications (with
[`ggiraph::girafeOutput()`](https://davidgohel.github.io/ggiraph/reference/girafeOutput.html)
and
[`ggiraph::renderGirafe()`](https://davidgohel.github.io/ggiraph/reference/renderGirafe.html)),
and
[`jp_save()`](https://skhademioore.github.io/joinpoint/reference/jp_save.md)
writes them to a standalone web page (without Pandoc, the page needs the
folder of scripts written beside it):

``` r

jp_save(jp_interactive(fit), file.path(dir, "mortality.html"))
```

By default the text of the graphic uses fonts found on the computer that
opens it, which keeps the page small; `embed_fonts = TRUE` embeds the
fonts, so that the graphic looks exactly the same everywhere, at the
cost of several megabytes.
