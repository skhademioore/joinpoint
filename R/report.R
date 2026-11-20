# Export of tables, figures and reports to document formats: delimited
# text, Markdown and LaTeX are written directly; Word, OpenDocument,
# PowerPoint, RTF, HTML and PDF documents are made by Pandoc (through the
# rmarkdown package).

#' Save a table, figure or report
#'
#' Writes a result table, a figure or an interactive graphic to a file whose
#' format follows its extension, for a manuscript, a report or a
#' presentation.
#'
#' * Tables (data frames such as the results of [jp_table()], [apc()] or
#'   [model_selection()]`$criteria`): `.csv`, `.tsv`, `.txt`, `.md`
#'   (Markdown pipe table), `.tex` (LaTeX `table` with booktabs rules; load
#'   `\usepackage{booktabs}`), `.html`, `.docx` (Word table), `.odt`,
#'   `.rtf`, `.pptx` and `.pdf`.
#' * Figures (ggplot objects such as [autoplot()][autoplot.joinpoint] or
#'   [plot_apc()]): the image formats of [ggplot2::ggsave()] (`.png`, `.pdf`,
#'   `.svg`, `.jpeg`, `.tiff`, `.eps`, ...), `.tex` (a LaTeX `figure` that
#'   includes a PDF of the figure written beside it), and `.html`, `.docx`,
#'   `.odt`, `.rtf`, `.pptx` and `.md` documents holding the figure.
#' * Interactive graphics ([jp_interactive()]): `.html`, a standalone web
#'   page.
#' * A `joinpoint` fit: a complete report, see [jp_report()].
#'
#' Word, OpenDocument, PowerPoint, RTF, HTML and PDF documents are made with
#' Pandoc, which comes with RStudio and Quarto (or see
#' <https://pandoc.org>), through the rmarkdown package; PDF documents also
#' need LaTeX (for example `tinytex::install_tinytex()`).
#'
#' @param x A data frame, a ggplot, an htmlwidget or a `joinpoint` fit.
#' @param file Output file; the extension sets the format.
#' @param title Optional caption of the table or figure.
#' @param note Optional note printed under the table or figure.
#' @param digits Decimals of the numeric columns of a table (character
#'   columns, such as those of [jp_table()], are written as they are).
#' @param width,height,dpi Size of a figure in inches and its resolution.
#' @param label Optional LaTeX label (`.tex` output).
#' @param ... Further arguments to [jp_report()] when `x` is a fit.
#' @return The path of the written file, invisibly.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' dir <- file.path(tempdir(), "joinpoint-save")
#' dir.create(dir)
#' jp_save(jp_table(fit), file.path(dir, "trends.tex"),
#'         title = "Trends in tuberculosis mortality")
#' jp_save(jp_table(fit), file.path(dir, "trends.md"))
#' jp_save(autoplot(fit), file.path(dir, "trend.png"), width = 6, height = 4)
#' if (requireNamespace("rmarkdown", quietly = TRUE) &&
#'     rmarkdown::pandoc_available()) {
#'   jp_save(jp_table(fit), file.path(dir, "trends.docx"))
#' }
#' list.files(dir)
#' unlink(dir, recursive = TRUE)
#' @export
jp_save <- function(x, file, title = NULL, note = NULL, digits = 3,
                    width = 7.5, height = 5, dpi = 300, label = NULL, ...) {
  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    stop("`file` must be a file name.", call. = FALSE)
  }
  ext <- jp_ext(file)
  if (inherits(x, "joinpoint")) {
    return(jp_report(x, file, title = title %||% "Joinpoint regression analysis",
                     ...))
  }
  if (inherits(x, "htmlwidget")) {
    if (ext != "html") {
      stop("Interactive graphics are saved as .html files; save the static ",
           "figure (interactive = FALSE) for other formats.", call. = FALSE)
    }
    jp_need("htmlwidgets", "saving interactive graphics")
    out <- jp_abs(file)
    htmlwidgets::saveWidget(x, out, title = title %||% "joinpoint",
                            selfcontained = jp_pandoc_ok(FALSE))
    return(invisible(file))
  }
  if (inherits(x, "ggplot")) {
    return(invisible(jp_save_figure(x, file, ext, title, note, width, height,
                                    dpi, label)))
  }
  if (is.data.frame(x)) {
    return(invisible(jp_save_table(x, file, ext, title, note, digits, label)))
  }
  stop("`x` must be a data frame, a ggplot, an htmlwidget or a joinpoint ",
       "fit.", call. = FALSE)
}

# Absolute path of a file that need not exist yet (in an existing folder).
jp_abs <- function(file) {
  file.path(normalizePath(dirname(file), mustWork = TRUE), basename(file))
}

jp_ext <- function(file) {
  ext <- tolower(sub(".*\\.", "", basename(file)))
  if (!grepl(".", basename(file), fixed = TRUE) || !nzchar(ext)) {
    stop("`file` needs an extension that names the format.", call. = FALSE)
  }
  ext
}

jp_doc_formats <- c("html", "docx", "odt", "rtf", "pptx", "pdf")

jp_save_table <- function(x, file, ext, title, note, digits, label) {
  tab <- jp_format_table(x, digits)
  switch(ext,
    csv = utils::write.csv(tab, file, row.names = FALSE),
    tsv = ,
    txt = utils::write.table(tab, file, sep = "\t", quote = FALSE,
                             row.names = FALSE),
    tex = writeLines(jp_latex_table(tab, title, note, label,
                                    jp_numeric_cols(x)), file),
    md = writeLines(jp_md_table(tab, title, note, jp_numeric_cols(x)), file),
    if (ext %in% jp_doc_formats) {
      jp_pandoc(jp_md_table(tab, title, note, jp_numeric_cols(x)), file, ext)
    } else {
      stop("Tables can be saved as .csv, .tsv, .txt, .md, .tex, .html, ",
           ".docx, .odt, .rtf, .pptx or .pdf files.", call. = FALSE)
    }
  )
  file
}

jp_save_figure <- function(x, file, ext, title, note, width, height, dpi,
                           label) {
  img <- c("png", "pdf", "svg", "jpeg", "jpg", "tiff", "tif", "bmp", "eps",
           "ps")
  if (ext %in% img) {
    ggplot2::ggsave(file, x, width = width, height = height, dpi = dpi,
                    bg = "white")
  } else if (ext == "tex") {
    pdf <- sub("\\.[^.]*$", ".pdf", file)
    ggplot2::ggsave(pdf, x, width = width, height = height, bg = "white")
    writeLines(c(
      "\\begin{figure}[htbp]", "\\centering",
      paste0("\\includegraphics[width=\\linewidth]{",
             sub("\\.pdf$", "", basename(pdf)), "}"),
      if (!is.null(title)) paste0("\\caption{", jp_tex(title), "}"),
      if (!is.null(label)) paste0("\\label{", label, "}"),
      if (!is.null(note)) paste0("\\par\\smallskip{\\footnotesize ",
                                 jp_tex(note), "}"),
      "\\end{figure}"
    ), file)
  } else if (ext %in% c(jp_doc_formats, "md")) {
    # The image goes beside a Markdown file, inside other documents.
    png <- if (ext == "md") {
      sub("\\.[^.]*$", ".png", file)
    } else {
      tempfile("figure", fileext = ".png")
    }
    ggplot2::ggsave(png, x, width = width, height = height, dpi = dpi,
                    bg = "white")
    md <- c(paste0("![", jp_md_escape(title %||% ""), "](",
                   if (ext == "md") basename(png) else png, ")"),
            if (!is.null(note)) c("", jp_md_escape(note)))
    if (ext == "md") writeLines(md, file) else jp_pandoc(md, file, ext)
  } else {
    stop("Figures can be saved as image files (.png, .pdf, .svg, ...), ",
         ".tex, .md, .html, .docx, .odt, .rtf or .pptx files.", call. = FALSE)
  }
  file
}

# Numeric columns of a table, formatted with fixed decimals (p values as
# "< 0.001"); missing values become empty cells.
jp_format_table <- function(x, digits) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  rownames(x) <- NULL
  for (j in seq_along(x)) {
    v <- x[[j]]
    if (is.numeric(v)) {
      pv <- grepl("^p(_value|_unadjusted)?$", names(x)[j])
      out <- if (pv) {
        ifelse(v < 0.001, "< 0.001", formatC(v, digits = 3, format = "f"))
      } else if (all(is.na(v) | v == round(v))) {
        format(v, scientific = FALSE, trim = TRUE)
      } else {
        formatC(v, digits = digits, format = "f")
      }
      out[is.na(v)] <- ""
      x[[j]] <- out
    } else {
      out <- as.character(v)
      out[is.na(out)] <- ""
      x[[j]] <- out
    }
  }
  x
}

jp_numeric_cols <- function(x) {
  vapply(x, is.numeric, logical(1)) |
    vapply(x, function(v) {
      v <- as.character(v[!is.na(v)])
      length(v) > 0L && all(grepl("^[-<>0-9. ,()*%]+$", v))
    }, logical(1))
}

jp_tex <- function(v) {
  v <- as.character(v)
  v <- gsub("\\", "\\textbackslash{}", v, fixed = TRUE)
  v <- gsub("([&%$#_{}])", "\\\\\\1", v)
  v <- gsub("~", "\\textasciitilde{}", v, fixed = TRUE)
  v <- gsub("^", "\\textasciicircum{}", v, fixed = TRUE)
  v <- gsub("<", "$<$", v, fixed = TRUE)
  gsub(">", "$>$", v, fixed = TRUE)
}

jp_latex_table <- function(tab, title = NULL, note = NULL, label = NULL,
                           right = rep(FALSE, ncol(tab))) {
  row <- function(v) paste0(paste(jp_tex(v), collapse = " & "), " \\\\")
  body <- vapply(seq_len(nrow(tab)), function(i) {
    row(unlist(tab[i, ], use.names = FALSE))
  }, character(1))
  c("\\begin{table}[htbp]", "\\centering",
    if (!is.null(title)) paste0("\\caption{", jp_tex(title), "}"),
    if (!is.null(label)) paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", paste(ifelse(right, "r", "l"), collapse = ""),
           "}"),
    "\\toprule", row(names(tab)), "\\midrule", body, "\\bottomrule",
    "\\end{tabular}",
    if (!is.null(note)) paste0("\\par\\smallskip{\\footnotesize ",
                               jp_tex(note), "}"),
    "\\end{table}")
}

jp_md_escape <- function(v) {
  v <- gsub("\\", "\\\\", as.character(v), fixed = TRUE)
  v <- gsub("|", "\\|", v, fixed = TRUE)
  v <- gsub("*", "\\*", v, fixed = TRUE)
  v <- gsub("_", "\\_", v, fixed = TRUE)
  gsub("<", "\\<", v, fixed = TRUE)
}

jp_md_table <- function(tab, title = NULL, note = NULL,
                        right = rep(FALSE, ncol(tab))) {
  row <- function(v) paste0("| ", paste(jp_md_escape(v), collapse = " | "),
                            " |")
  c(row(names(tab)),
    paste0("|", paste(ifelse(right, "---:", ":---"), collapse = "|"), "|"),
    vapply(seq_len(nrow(tab)), function(i) {
      row(unlist(tab[i, ], use.names = FALSE))
    }, character(1)),
    if (!is.null(title)) c("", paste0("Table: ", jp_md_escape(title))),
    if (!is.null(note)) c("", jp_md_escape(note)))
}

jp_pandoc_ok <- function(error = TRUE) {
  ok <- requireNamespace("rmarkdown", quietly = TRUE) &&
    rmarkdown::pandoc_available()
  if (!ok && error) {
    stop("Pandoc is needed for this format; it comes with RStudio and ",
         "Quarto, or see https://pandoc.org (and install the rmarkdown ",
         "package).", call. = FALSE)
  }
  ok
}

# Converts Markdown lines to the document `file` with Pandoc.
jp_pandoc <- function(md, file, ext) {
  jp_pandoc_ok()
  src <- tempfile("joinpoint", fileext = ".md")
  on.exit(unlink(src))
  writeLines(md, src, useBytes = TRUE)
  opts <- c("--standalone", "--metadata", "pagetitle=joinpoint")
  if (ext == "html") {
    opts <- c(opts, if (rmarkdown::pandoc_available("2.19")) {
      "--embed-resources"
    } else {
      "--self-contained"
    })
  }
  out <- jp_abs(file)
  res <- tryCatch(
    rmarkdown::pandoc_convert(src, output = out, from = "markdown",
                              options = opts),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    stop("Pandoc could not write ", basename(file),
         if (ext == "pdf") " (PDF output needs LaTeX, for example tinytex)",
         ": ", conditionMessage(res), call. = FALSE)
  }
  invisible(file)
}

#' Joinpoint report as a document, R Markdown or Quarto file
#'
#' Writes a report of a joinpoint fit: the settings of the analysis, the
#' trends of the selected models (APCs and AAPCs with confidence intervals),
#' the joinpoints, the model selection statistics and the figures. The
#' format follows the extension of `file`:
#'
#' * `.Rmd` or `.qmd`: an R Markdown or Quarto source document, with the fit
#'   saved beside it (`<name>_fit.rds`) and read by its first chunk, ready to
#'   edit and render (in RStudio, with [rmarkdown::render()], or with
#'   `quarto render`).
#' * `.html`, `.docx`, `.odt`, `.pdf`, `.md` or `.tex`: the rendered
#'   report (with the rmarkdown package and Pandoc; PDF also needs LaTeX).
#'   The figures of HTML reports are interactive when the ggiraph package is
#'   installed; the LaTeX and Markdown reports keep their figures in a
#'   folder beside them.
#'
#' @param object A `joinpoint` object.
#' @param file Output file; the extension sets the format.
#' @param title,author Title and author of the report.
#' @param model `"final"` or a number of joinpoints.
#' @param plots Figures to include: any of `"trend"`, `"apc"`, `"aapc"` and
#'   `"selection"`, or none (`character()`).
#' @param format Output format written into the header of an `.Rmd` or
#'   `.qmd` source: `"html"`, `"docx"` or `"pdf"`.
#' @param interactive Interactive figures (HTML output only); by default
#'   when the output is HTML and the ggiraph package is installed.
#' @param quiet Render without printing the progress of knitr and Pandoc.
#' @return The path of the written file, invisibly.
#' @examples
#' tb <- subset(us_mortality, cause == "Tuberculosis")
#' fit <- joinpoint(tb, year, rate, kmax = 2, ci_method = "parametric")
#' dir <- file.path(tempdir(), "joinpoint-report")
#' dir.create(dir)
#' jp_report(fit, file.path(dir, "tuberculosis.qmd"))
#' list.files(dir)
#' \donttest{
#' if (requireNamespace("rmarkdown", quietly = TRUE) &&
#'     rmarkdown::pandoc_available()) {
#'   jp_report(fit, file.path(dir, "tuberculosis.docx"))
#' }
#' }
#' unlink(dir, recursive = TRUE)
#' @export
jp_report <- function(object, file, title = "Joinpoint regression analysis",
                      author = NULL, model = "final",
                      plots = c("trend", "apc", "aapc", "selection"),
                      format = c("html", "docx", "pdf"), interactive = NULL,
                      quiet = TRUE) {
  stopifnot(inherits(object, "joinpoint"))
  ext <- jp_ext(file)
  format <- match.arg(format)
  if (length(plots)) {
    plots <- match.arg(plots, c("trend", "apc", "aapc", "selection"),
                       several.ok = TRUE)
  }
  if (is.null(aapc(object, model))) plots <- setdiff(plots, "aapc")
  out_format <- switch(ext, rmd = format, qmd = format, ext)
  if (!out_format %in% c("html", "docx", "odt", "pdf", "md", "tex")) {
    stop("Reports are written as .Rmd, .qmd, .html, .docx, .odt, .pdf, .md ",
         "or .tex files.", call. = FALSE)
  }
  interactive <- interactive %||%
    (out_format == "html" && requireNamespace("ggiraph", quietly = TRUE))
  if (interactive && out_format != "html") {
    stop("Interactive figures need HTML output.", call. = FALSE)
  }
  if (ext %in% c("rmd", "qmd")) {
    rds <- sub("\\.[^.]*$", "_fit.rds", file)
    saveRDS(object, rds)
    writeLines(jp_report_source(object, basename(rds), title, author, model,
                                plots, out_format, interactive,
                                quarto = ext == "qmd"), file)
    return(invisible(file))
  }
  jp_pandoc_ok()
  wd <- tempfile("joinpoint-report")
  dir.create(wd)
  on.exit(unlink(wd, recursive = TRUE))
  saveRDS(object, file.path(wd, "fit.rds"))
  src <- file.path(wd, "report.Rmd")
  writeLines(jp_report_source(object, "fit.rds", title, author, model, plots,
                              out_format, interactive, quarto = FALSE), src)
  of <- switch(out_format,
    html = "rmarkdown::html_document(self_contained = TRUE)",
    docx = "rmarkdown::word_document()",
    odt = "rmarkdown::odt_document()",
    pdf = "rmarkdown::pdf_document()",
    md = "rmarkdown::md_document(variant = \"gfm\")",
    tex = "rmarkdown::latex_document()"
  )
  out <- jp_abs(file)
  # Rendered by a new R session, which leaves the current session (and a
  # document being knitted, when the report is made from one) untouched.
  q <- function(v) deparse(v, width.cutoff = 500L)
  script <- file.path(wd, "render.R")
  writeLines(c(
    paste0(".libPaths(", q(.libPaths()), ")"),
    paste0("rmarkdown::find_pandoc(dir = ", q(rmarkdown::find_pandoc()$dir),
           ")"),
    paste0("rmarkdown::render(", q(src), ", output_format = ", of,
           ", output_file = ", q(basename(out)), ", output_dir = ",
           q(dirname(out)), ", intermediates_dir = ", q(wd),
           ", knit_root_dir = ", q(wd), ", quiet = ", q(isTRUE(quiet)), ")")
  ), script)
  log <- file.path(wd, "render.log")
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--no-save", "--no-restore", shQuote(script)),
                    stdout = if (isTRUE(quiet)) log else "",
                    stderr = if (isTRUE(quiet)) log else "")
  if (!identical(as.integer(status), 0L) || !file.exists(out)) {
    msg <- if (file.exists(log)) utils::tail(readLines(log), 12L) else ""
    stop("The report could not be rendered",
         if (out_format == "pdf") " (PDF output needs LaTeX, for example tinytex)",
         ":\n", paste(msg, collapse = "\n"), call. = FALSE)
  }
  invisible(file)
}

# Source of the report: a YAML header, a setup chunk reading the fit, the
# settings as text, the tables and the figures.
jp_report_source <- function(object, rds, title, author, model, plots,
                             format, interactive, quarto) {
  s <- object$settings
  nm <- jp_stat_names(s$interval, s$log, s$family)
  q <- function(v) paste0('"', gsub('"', '\\\\"', v), '"')
  yaml <- if (quarto) {
    c("---", paste("title:", q(title)),
      if (!is.null(author)) paste("author:", q(author)),
      paste("date:", q(format(Sys.Date()))),
      paste0("format: ", switch(format, html = "\n  html:\n    embed-resources: true",
                                format)),
      "execute:", "  echo: false", "  warning: false", "  message: false",
      "---")
  } else {
    c("---", paste("title:", q(title)),
      if (!is.null(author)) paste("author:", q(author)),
      paste("date:", q(format(Sys.Date()))),
      paste("output:", switch(format, html = "html_document",
                              docx = "word_document", pdf = "pdf_document",
                              odt = "odt_document", md = "md_document",
                              tex = "latex_document")),
      "---")
  }
  n_coh <- sum(!vapply(object$fits, is.null, logical(1)))
  rows <- ceiling(n_coh / ceiling(sqrt(n_coh)))
  fh <- round(min(9, 1.2 + 2.6 * rows), 1)
  fig <- function(code, height = fh) {
    c(paste0("```{r, fig.width = 7.5, fig.height = ", height, "}"), code,
      "```", "")
  }
  plot_code <- function(p) {
    fun <- switch(p, trend = "autoplot", apc = "plot_apc", aapc = "plot_aapc",
                  selection = "plot_selection")
    margs <- if (p %in% c("trend", "apc", "aapc") && !identical(model, "final")) {
      paste0(", model = ", model)
    } else {
      ""
    }
    if (interactive) {
      paste0("jp_interactive(fit, plot = \"", p, "\"", margs, ")")
    } else {
      paste0(fun, "(fit", margs, ")")
    }
  }
  mstr <- if (identical(model, "final")) "" else paste0(", model = ", model)
  ks <- vapply(Filter(Negate(is.null), object$fits), `[[`, numeric(1),
               "k_selected")
  method <- switch(s$method, wbic = "the weighted BIC", bic = "the BIC",
                   bic3 = "BIC3", mbic = "the modified BIC",
                   dds = "data dependent selection",
                   permutation = "the permutation test",
                   lrt = "the Monte Carlo likelihood ratio test")
  model_txt <- switch(s$family,
    gaussian = paste0("a ", if (s$log) "log-linear" else "linear",
                      " joinpoint model with ", s$variance, " variance"),
    poisson = "a Poisson joinpoint regression with log link",
    binomial = "a logistic joinpoint regression"
  )
  txt <- paste0(
    "The trend of ", jp_md_escape(jp_title(s$y %||% "the response")),
    " over ", jp_md_escape(jp_title(s$x)), " was modeled by ", model_txt,
    if (n_coh > 1L) paste0(" for each of ", n_coh, " cohorts") else "",
    ". The number of joinpoints (", s$kmin, " to at most ",
    max(vapply(Filter(Negate(is.null), object$fits), `[[`, numeric(1),
               "kmax")), ") was selected by ", method, ", and ",
    if (s$ci_method == "empirical") "empirical quantile" else "parametric",
    " confidence intervals were computed. The selected models have ",
    paste(sort(unique(ks)), collapse = ", "), " joinpoint",
    if (all(ks == 1)) "" else "s", "."
  )
  c(yaml, "",
    "```{r setup, include = FALSE}",
    if (!quarto) paste0("knitr::opts_chunk$set(echo = FALSE, warning = FALSE, ",
                        "message = FALSE, dpi = 150)"),
    "library(joinpoint)",
    paste0("fit <- readRDS(\"", rds, "\")"),
    "```", "",
    "## Analysis", "", txt, "",
    paste0("## ", nm$apc, "s"), "",
    "```{r}",
    paste0("knitr::kable(jp_table(fit", mstr, "), caption = \"",
           nm$apc_abbr, "s of the segments and ", nm$aapc_abbr,
           "s, with 95% confidence intervals.\")"),
    "```", "",
    if (any(ks > 0)) c(
      "## Joinpoints", "",
      "```{r}",
      paste0("knitr::kable(subset(joinpoints(fit", mstr, "), select = ",
             "-c(model, final, ci_method)), digits = 2, ",
             "caption = \"Estimated joinpoints with 95% confidence ",
             "intervals.\")"),
      "```", ""),
    "## Model selection", "",
    "```{r}",
    paste0("knitr::kable(model_selection(fit)$criteria, digits = 3, ",
           "caption = \"Selection statistics of the fitted models.\")"),
    "```", "",
    if (length(plots)) c("## Figures", "",
                         unlist(lapply(plots, function(p) {
                           fig(plot_code(p), if (p == "selection") {
                             min(fh, 4.5)
                           } else if (p == "aapc") {
                             round(min(9, 2 + 0.35 * n_coh), 1)
                           } else {
                             fh
                           })
                         }))))
}

#' Printing joinpoint fits in R Markdown and Quarto
#'
#' A fit printed in a chunk of an R Markdown or Quarto document shows the
#' same summary as in the console. With `options(joinpoint.knit_table =
#' TRUE)` it appears instead as the publication table of [jp_table()],
#' formatted for the output (HTML, LaTeX, Word or Markdown); the chunk option
#' `tab.cap` sets its caption.
#'
#' @param x A `joinpoint` object.
#' @param options Chunk options (passed by knitr).
#' @param ... Further arguments to the table method.
#' @return The output of [knitr::knit_print()].
#' @keywords internal
#' @exportS3Method knitr::knit_print
knit_print.joinpoint <- function(x, options = NULL, ...) {
  if (!isTRUE(getOption("joinpoint.knit_table", FALSE))) {
    return(knitr::normal_print(x))
  }
  nm <- jp_stat_names(x$settings$interval, x$settings$log, x$settings$family)
  cap <- options$tab.cap %||% paste0(nm$apc_abbr, "s and ", nm$aapc_abbr,
                                     "s with 95% confidence intervals.")
  knitr::knit_print(knitr::kable(jp_table(x), caption = cap), options = options,
                    ...)
}
