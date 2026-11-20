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

out_fit <- function(piece = TRUE, empirical = FALSE, custom = FALSE) {
  x <- 1:20
  d <- data.frame(x, y = exp(2 + .04 * x + if (piece) {
    -.2 * pmax(x - 10, 0) + .015 * sin(x)
  } else {
    .1 * sin(x)
  }))
  ctl <- if (custom) {
    jp_control(alpha_apc = .1, alpha_aapc = .2, alpha_joinpoint = .01)
  } else {
    jp_control(n_resample = 99)
  }
  joinpoint(d, x, y, kmax = 1, method = "bic",
            ci_method = if (empirical) "empirical" else "parametric",
            control = ctl)
}

test_that("outputs state the configured confidence levels", {
  f <- out_fit(custom = TRUE)
  dir <- temp_dir()
  jp_report(f, file.path(dir, "levels.Rmd"), interactive = FALSE)
  src <- readLines(file.path(dir, "levels.Rmd"))
  expect_false(any(grepl("95% confidence", src, fixed = TRUE)))
  expect_true(any(grepl("90% confidence intervals for the APCs and 80%", src,
                        fixed = TRUE)))
  expect_true(any(grepl("99% confidence intervals", src, fixed = TRUE)))
  skip_if_not_installed("ggiraph")
  b <- ggplot2::ggplot_build(autoplot(f, interactive = TRUE))
  tips <- unlist(lapply(b$data, function(d) d$tooltip), use.names = FALSE)
  sg <- tips[grepl("Segment ", tips, fixed = TRUE)]
  expect_true(length(sg) > 0 && all(grepl("90% CI", sg, fixed = TRUE)))
  jt <- tips[grepl("Joinpoint ", tips, fixed = TRUE)]
  expect_true(all(grepl("99% CI", jt, fixed = TRUE)))
})

test_that("reports describe the requested model", {
  f <- out_fit(piece = FALSE, empirical = TRUE)
  expect_equal(f$fits[[1]]$k_selected, 0)
  dir <- temp_dir()
  jp_report(f, file.path(dir, "m1.Rmd"), model = 1, plots = character(),
            interactive = FALSE)
  src <- readLines(file.path(dir, "m1.Rmd"))
  expect_true("## Joinpoints" %in% src)
  expect_false(any(grepl("empirical quantile", src, fixed = TRUE)))
  expect_true(any(grepl("models with 1 joinpoint, as requested", src,
                        fixed = TRUE)))
  expect_error(jp_report(f, file.path(dir, "m3.Rmd"), model = 3),
               "No fitted model")
  expect_error(jp_report(f, file.path(dir, "m.Rmd"), model = "x"), "model")
})

test_that("a report of a model without joinpoints renders", {
  skip_on_cran()
  skip_if_not(has_pandoc())
  f <- out_fit()
  expect_equal(f$fits[[1]]$k_selected, 1)
  dir <- temp_dir()
  jp_report(f, file.path(dir, "m0.md"), model = 0, plots = character(),
            interactive = FALSE)
  md <- readLines(file.path(dir, "m0.md"))
  expect_false("## Joinpoints" %in% md)
})

test_that("titles and file names are written as valid YAML and R", {
  f <- out_fit()
  dir <- temp_dir()
  title <- "Trend C:\\rates \"quoted\"\nsecond line"
  jp_report(f, file.path(dir, "meta.Rmd"), title = title, plots = character(),
            interactive = FALSE)
  if (requireNamespace("rmarkdown", quietly = TRUE)) {
    expect_equal(rmarkdown::yaml_front_matter(file.path(dir, "meta.Rmd"))$title,
                 title)
  }
  path <- file.path(dir, "quoted \"study\".Rmd")
  jp_report(f, path, plots = character(), interactive = FALSE)
  code <- grep("readRDS", readLines(path), value = TRUE, fixed = TRUE)
  expect_equal(eval(parse(text = sub("fit <- readRDS", "identity", code))),
               "quoted \"study\"_fit.rds")
})

test_that("Markdown reports link their figures by relative paths", {
  skip_on_cran()
  skip_if_not(has_pandoc())
  f <- out_fit(piece = FALSE)
  dir <- temp_dir()
  jp_report(f, file.path(dir, "rep.md"), plots = "trend", interactive = FALSE)
  src <- readLines(file.path(dir, "rep.md"))
  refs <- sub(".*!\\[[^]]*\\]\\(([^)]*)\\).*", "\\1", grep("!\\[", src, value = TRUE))
  expect_gt(length(refs), 0)
  expect_false(any(grepl("^/|^[A-Za-z]:|^file:", refs)))
  expect_true(all(file.exists(file.path(dir, refs))))
})

test_that("cell contents survive tab separated and LaTeX tables", {
  dir <- temp_dir()
  tab <- data.frame(group = c("A\tB", "C\nD"), estimate = c(1, 2))
  jp_save(tab, file.path(dir, "d.tsv"))
  back <- utils::read.delim(file.path(dir, "d.tsv"), check.names = FALSE)
  expect_equal(back$group, tab$group)
  expect_equal(joinpoint:::jp_tex("a\\b {c} ~"),
               "a\\textbackslash{}b \\{c\\} \\textasciitilde{}")
})
