L <- function(d, ...) {
  joinpoint(d, x, y, log = FALSE, type = "other", ci_method = "parametric",
            ...)
}

test_that("cohort keys are matched as tuples", {
  d <- data.frame(x = 1:10, y = 1:10 + sin(1:10),
                  g = rep(c("A.B", "A"), each = 5),
                  h = rep(c("C", "B.C"), each = 5))
  f <- L(d, by = c("g", "h"), kmax = 0)
  expect_equal(nrow(f$cohorts), 2L)
  expect_equal(vapply(f$fits, function(g) length(g$x), integer(1)), c(5L, 5L),
               ignore_attr = TRUE)
  a <- data.frame(g = c("a.1", "a"), time = c(2, 1.2), age = "all",
                  count = c(10, 90), pop = 100)
  out <- age_adjust(a, time, count, pop, age, c(all = 1), by = "g",
                    rate_per = 100)
  expect_equal(out$rate[match(a$g, out$g)], c(10, 90))
})

test_that("by-variables cannot take the names of result columns", {
  d <- data.frame(x = 1:6, y = 1:6 + sin(1:6), fit = "a")
  expect_error(L(d, by = "fit", kmax = 0), "rename fit")
})

test_that("comparisons validate their arguments and observation times", {
  d <- data.frame(x = rep(1:5, 2), y = rep(1:5, 2) + rep(c(0, .1), each = 5) *
                    sin(1:5), g = rep(c("a", "b"), each = 5))
  f <- L(d, by = "g", kmax = 0)
  expect_error(jp_compare(f, n_perm = 1.5), "n_perm")
  expect_error(jp_compare(f, alpha = 2), "alpha")
  d2 <- data.frame(x = c(1e9 + 1:5, 1e9 + 2:6), y = d$y, g = d$g)
  expect_error(jp_compare(L(d2, by = "g", kmax = 0), n_perm = 9),
               "same independent")
})

test_that("predictions keep the rows of newdata", {
  d <- data.frame(x = rep(1:5, 2), y = rep(1:5, 2) + rep(c(0, 1), each = 5),
                  g = rep(c("a", "b"), each = 5))
  f <- L(d, by = "g", kmax = 0)
  p <- predict(f, data.frame(x = c(1, 2, 3), g = c("a", NA, "c")))
  expect_equal(nrow(p), 3L)
  expect_equal(p$fit, c(1, NA, NA))
  nd <- data.frame(x = c(4, 2, 3), g = c("b", "a", "b"))
  p <- predict(f, nd)
  expect_equal(p$g, nd$g)
  expect_equal(p$x, nd$x)
  expect_equal(p$fit, c(5, 2, 4))
})

one_k <- function(f) as.character(f$fits[[1]]$k_selected)

test_that("the all-model table identifies each model", {
  f <- L(data.frame(x = 1:20, y = sin(1:20)), kmax = 1)
  tab <- jp_table(f, model = "all")
  expect_true(all(c("Model", "Final") %in% names(tab)))
  expect_equal(sum(tab$Trend == "AAC"), 2L)
  expect_equal(tab$Final[tab$Model == one_k(f)][1], "Yes")
})

test_that("non-significant results are not labeled as equality", {
  d <- data.frame(x = rep(1:12, 2), g = rep(c("a", "b"), each = 12), y = 1)
  f <- L(d, by = "g", kmax = 0)
  expect_false(any(c("Stable", "Parallel") %in% jp_compare(f, n_perm = 19)$result))
  expect_equal(as.character(joinpoint:::jp_trend(0, 1, .05)),
               "No significant change")
})

test_that("binomial jumps are odds ratios", {
  d <- data.frame(x = 1:20, n = rep(c(2000, 4000), each = 10), p = 10000)
  f <- joinpoint(d, x, count = n, population = p, family = "binomial",
                 jump = 10, kmax = 0, method = "bic")
  j <- jump_estimates(f)
  expect_equal(j$estimate[2], 8 / 3, tolerance = 1e-7)
  expect_equal(j$statistic[2], "Odds ratio")
})

test_that("exports replace the files of an earlier export", {
  x <- 1:12
  d <- data.frame(x, y = 10 + x + pmax(x - 6, 0) + .1 * sin(x))
  dir <- tempfile()
  on.exit(unlink(dir, recursive = TRUE))
  jp_export(joinpoint(d, x, y, kmin = 1, kmax = 1, ci_method = "parametric"),
            dir)
  jp_export(joinpoint(d, x, y, kmax = 0, ci_method = "parametric"), dir)
  expect_false(file.exists(file.path(dir, "joinpoint.Joinpoints.txt")))
})

test_that("cluster rates standardize the original or compressed age groups", {
  d <- expand.grid(x = 1:4, a = 1:2)
  d$p <- ifelse(d$a == 1, 100 * d$x, 100)
  d$n <- ifelse(d$a == 1, 0, 10)
  std <- c("1" = 1, "2" = 1)
  cl <- jp_cluster(d, x, n, p, a, std, kmax = 0, rate_per = 1,
                   ci_method = "parametric")
  expect_equal(cl$rates$rate, age_adjust(d, x, n, p, a, std,
                                         rate_per = 1)$rate)
  # Compressed: the pooled rate 10 / (100 x + 100).
  cc <- jp_cluster(d, x, n, p, a, std, kmax = 0, rate_per = 1,
                   standardize = "compressed", ci_method = "parametric")
  expect_equal(cc$rates$rate, 10 / (100 * (1:4) + 100))
})

test_that("clustering validates its controls and bounds its search", {
  d <- expand.grid(x = 1:20, a = 1:3)
  d$p <- 1e4
  d$n <- 100 + d$x
  std <- c("1" = 1, "2" = 1, "3" = 1)
  for (m in list(0, -1, 1.5, NA, c(1, 2))) {
    expect_error(jp_cluster(d, x, n, p, a, std, max_clusters = m),
                 "max_clusters")
  }
  suppressMessages(expect_message(
    jp_cluster(d, x, n, p, a, std, ci_method = "parametric", kmin = 0,
               control = jp_control(time_budget = 1e-6)),
    "Clustering: recommended maximum"
  ))
})

test_that("cohorts with the same display label stay apart", {
  mk <- function(sep) {
    d <- expand.grid(x = 1:8, g = 1:2)
    d$a <- ifelse(d$g == 1, paste0("A", sep, "B"), "A")
    d$b <- ifelse(d$g == 1, "C", paste0("B", sep, "C"))
    d$y <- exp(1 + .05 * d$x + log(10) * (d$g == 2) + .01 * sin(d$x))
    d
  }
  d <- mk(", ")
  f <- joinpoint(d, x, y, by = c("a", "b"), kmax = 0, ci_method = "parametric")
  b <- ggplot2::ggplot_build(autoplot(f, labels = "none"))
  expect_equal(nrow(b$layout$layout), 2L)
  expect_equal(length(unique(plot_aapc(f)$data$cohort)), 2L)
  # Prediction matches tuples, whatever characters the keys hold.
  d <- mk("\r")
  f <- joinpoint(d, x, y, by = c("a", "b"), kmax = 0, ci_method = "parametric")
  one <- vapply(1:2, function(g) {
    predict(joinpoint(d[d$g == g, ], x, y, kmax = 0, ci_method = "parametric"),
            data.frame(x = 5))$fit
  }, numeric(1))
  expect_equal(predict(f, d[c(5, 13), c("x", "a", "b")])$fit, one,
               tolerance = 1e-10)
})

test_that("a by-variable named Model is refused", {
  d <- data.frame(x = rep(1:12, 2), Model = rep(c("a", "b"), each = 12))
  d$y <- exp(2 + .03 * d$x + .01 * sin(d$x) + .3 * (d$Model == "b"))
  expect_error(joinpoint(d, x, y, by = "Model", kmax = 0), "rename Model")
})

test_that("exactly parallel age trends form one cluster", {
  d <- expand.grid(x = 1:20, age = 1:3)
  d$count <- 100 * d$age
  d$pop <- 10000
  cl <- jp_cluster(d, x, count, pop, age, c("1" = 1, "2" = 1, "3" = 1),
                   kmax = 0, ci_method = "parametric")
  expect_equal(cl$selection$clusters[cl$selection$selected], 1)
  expect_equal(cl$settings$standardize, "original")
})

test_that("excluded rows are not validated", {
  d <- data.frame(x = 1:12, n = c(-1, 5:15), p = c(0, rep(100, 11)))
  f <- joinpoint(d, x, count = n, population = p, exclude = 1, kmax = 0,
                 ci_method = "parametric")
  expect_true(is.na(f$fits[[1]]$data$r[1]))
  expect_equal(nobs(f), 11)
})
