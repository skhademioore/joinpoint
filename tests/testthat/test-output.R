tb <- subset(us_mortality, cause %in% c("Stroke", "Tuberculosis"))
fit <- joinpoint(tb, year, rate, by = "cause", kmax = 2,
                 ci_method = "parametric")
has_pandoc <- function() {
  requireNamespace("rmarkdown", quietly = TRUE) &&
    rmarkdown::pandoc_available()
}

test_that("interactive graphics carry tooltips for segments and points", {
  skip_if_not_installed("ggiraph")
  skip_if_not_installed("gdtools")
  g <- autoplot(fit, interactive = TRUE, error_bars = TRUE)
  b <- ggplot2::ggplot_build(g)
  tips <- unlist(lapply(b$data, function(d) d$tooltip))
  expect_true(any(grepl("Segment 1: 1900-", tips, fixed = TRUE)))
  expect_true(any(grepl("Observed", tips, fixed = TRUE)))
  expect_true(any(grepl("Joinpoint", tips, fixed = TRUE)))
  # The curve and the label of a segment share its data id.
  ids <- unlist(lapply(b$data, function(d) d$data_id))
  expect_gt(sum(ids == "Stroke|segment|1"), 1L)
  for (p in c("trend", "apc", "aapc", "selection")) {
    w <- jp_interactive(fit, plot = p)
    expect_s3_class(w, "girafe")
    expect_false(any(grepl("liberation",
                           vapply(w$dependencies, `[[`, "", "name"))))
  }
  expect_s3_class(jp_interactive(g + ggplot2::labs(title = "Mortality")),
                  "girafe")
  expect_error(jp_interactive(1), "joinpoint fit or a ggplot")
  expect_error(jp_interactive(fit, zoom = 0), "zoom")
})

test_that("static figures are unchanged by the interactive option", {
  g <- autoplot(fit)
  expect_false(any(vapply(g$layers, function(l) {
    grepl("Interactive", class(l$geom)[1])
  }, logical(1))))
})

test_that("tables are written as text, Markdown and LaTeX", {
  dir <- temp_dir()
  tab <- jp_table(fit)
  f <- file.path(dir, c("t.csv", "t.tsv", "t.md", "t.tex"))
  for (x in f) jp_save(tab, x, title = "Trends_by cause", note = "* p < 0.05.")
  expect_equal(nrow(utils::read.csv(f[1])), nrow(tab))
  md <- readLines(f[3])
  expect_match(md[1], "^\\| cause \\| Trend")
  expect_true(any(md == "Table: Trends\\_by cause"))
  tex <- readLines(f[4])
  expect_true(all(c("\\toprule", "\\midrule", "\\bottomrule") %in% tex))
  expect_true("\\caption{Trends\\_by cause}" %in% tex)
  expect_true(any(grepl("$<$ 0.001", tex, fixed = TRUE)))
  # Numeric tables: fixed decimals, p values and missing cells.
  a <- data.frame(name = c("a_b", "c"), est = c(1.23456, NA),
                  p_value = c(0.0001, 0.2))
  jp_save(a, file.path(dir, "a.tex"), digits = 2)
  tex <- readLines(file.path(dir, "a.tex"))
  expect_true("a\\_b & 1.23 & $<$ 0.001 \\\\" %in% tex)
  expect_true("c &  & 0.200 \\\\" %in% tex)
  expect_error(jp_save(tab, file.path(dir, "t.xyz")), "Tables can be saved")
  expect_error(jp_save(tab, file.path(dir, "noext")), "extension")
  expect_error(jp_save(1:3, file.path(dir, "x.csv")), "data frame")
})

test_that("figures are written as images and LaTeX figures", {
  dir <- temp_dir()
  g <- autoplot(fit)
  jp_save(g, file.path(dir, "f.png"), width = 4, height = 3, dpi = 72)
  jp_save(g, file.path(dir, "f.tex"), title = "Trends", label = "fig:trend")
  expect_true(file.exists(file.path(dir, "f.png")))
  expect_true(file.exists(file.path(dir, "f.pdf")))
  tex <- readLines(file.path(dir, "f.tex"))
  expect_true("\\includegraphics[width=\\linewidth]{f}" %in% tex)
  expect_true("\\label{fig:trend}" %in% tex)
  jp_save(g, file.path(dir, "g.md"), title = "Trends", width = 4, height = 3,
          dpi = 72)
  expect_equal(readLines(file.path(dir, "g.md")), "![Trends](g.png)")
  expect_error(jp_save(g, file.path(dir, "f.xlsx")), "Figures can be saved")
})

test_that("documents are made with Pandoc", {
  skip_on_cran()
  skip_if_not(has_pandoc())
  dir <- temp_dir()
  for (e in c("docx", "html", "odt")) {
    jp_save(jp_table(fit), file.path(dir, paste0("t.", e)))
    expect_true(file.exists(file.path(dir, paste0("t.", e))))
  }
  jp_save(autoplot(fit), file.path(dir, "f.docx"), width = 4, height = 3,
          dpi = 72)
  expect_true(file.exists(file.path(dir, "f.docx")))
})

test_that("interactive graphics are saved as web pages only", {
  skip_on_cran()
  skip_if_not_installed("ggiraph")
  skip_if_not_installed("gdtools")
  dir <- temp_dir()
  w <- jp_interactive(fit)
  jp_save(w, file.path(dir, "w.html"))
  expect_true(file.exists(file.path(dir, "w.html")))
  expect_error(jp_save(w, file.path(dir, "w.png")), "html")
})

test_that("reports are written as R Markdown and Quarto sources", {
  dir <- temp_dir()
  jp_report(fit, file.path(dir, "r.Rmd"), title = "Mortality \"trends\"",
            author = "A. Author", format = "docx", interactive = FALSE)
  src <- readLines(file.path(dir, "r.Rmd"))
  expect_true("output: word_document" %in% src)
  expect_true("title: \"Mortality \\\"trends\\\"\"" %in% src)
  expect_true("fit <- readRDS(\"r_fit.rds\")" %in% src)
  expect_true(any(grepl("plot_apc(fit)", src, fixed = TRUE)))
  f2 <- readRDS(file.path(dir, "r_fit.rds"))
  expect_equal(apc(f2), apc(fit))
  jp_report(fit, file.path(dir, "q.qmd"), plots = "trend", model = 1)
  q <- readLines(file.path(dir, "q.qmd"))
  expect_true("format: " %in% q || any(grepl("^format:", q)))
  expect_true(any(grepl("model = 1", q, fixed = TRUE)))
  expect_false(any(grepl("plot_apc", q, fixed = TRUE)))
  expect_error(jp_report(fit, file.path(dir, "r.xlsx")), "Reports are")
  expect_error(jp_report(fit, file.path(dir, "r.docx"), interactive = TRUE),
               "HTML")
})

test_that("reports are rendered to documents", {
  skip_on_cran()
  skip_if_not(has_pandoc())
  dir <- temp_dir()
  jp_report(fit, file.path(dir, "r.docx"), plots = "trend")
  expect_true(file.exists(file.path(dir, "r.docx")))
  jp_save(fit, file.path(dir, "s.md"), plots = character())
  md <- readLines(file.path(dir, "s.md"))
  expect_true(any(grepl("Segment 1", md, fixed = TRUE)))
})

test_that("printed fits in knitr documents are tables", {
  skip_if_not_installed("knitr")
  expect_output(knitr::knit_print(fit), "Joinpoint")
  old <- options(joinpoint.knit_table = TRUE)
  on.exit(options(old))
  out <- knitr::knit_print(fit)
  expect_match(as.character(out), "Segment 1", fixed = TRUE)
})
