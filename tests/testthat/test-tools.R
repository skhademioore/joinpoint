test_that("age adjustment matches the direct standardization formula", {
  d <- data.frame(t = rep(1:2, each = 3), age = rep(c("a", "b", "c"), 2),
                  n = c(5, 10, 20, 6, 12, 25), p = rep(c(1e4, 2e4, 5e3), 2))
  std <- c(a = 1, b = 2, c = 1)
  out <- age_adjust(d, t, n, p, age, std)
  w <- std / sum(std)
  r1 <- 1e5 * sum(w * d$n[1:3] / d$p[1:3])
  se1 <- 1e5 * sqrt(sum(w^2 * d$n[1:3] / d$p[1:3]^2))
  expect_equal(out$rate[1], r1)
  expect_equal(out$se[1], se1)
})

test_that("pairwise comparison accepts parallel trends", {
  set.seed(4)
  x <- 1:25
  mu <- 2 + 0.05 * x - 0.08 * pmax(x - 12, 0)
  dd <- rbind(
    data.frame(g = "a", x = x, y = exp(mu + rnorm(25, sd = 0.03))),
    data.frame(g = "b", x = x, y = exp(mu + 0.5 + rnorm(25, sd = 0.03)))
  )
  fit <- joinpoint(dd, x, y, by = "g", kmax = 2, method = "bic",
                   ci_method = "parametric")
  cmp <- jp_compare(fit, n_perm = 199)
  expect_equal(cmp$result, "Parallel")
  cmp2 <- jp_compare(fit, type = "coincidence", n_perm = 199)
  expect_equal(cmp2$result, "Not coincident")
  expect_s3_class(autoplot(cmp2), "ggplot")
})

test_that("clustering recovers contiguous age groups with common trends", {
  set.seed(3)
  ages <- sprintf("%02d-%02d", seq(0, 70, 10), seq(9, 79, 10))
  sim <- expand.grid(year = 1990:2019, age = ages, stringsAsFactors = FALSE)
  ai <- match(sim$age, ages)
  lr <- log(5) + 0.35 * ai +
    ifelse(ai <= 4, 0.02 * (sim$year - 1990),
           0.01 * (sim$year - 1990) - 0.04 * pmax(sim$year - 2005, 0))
  sim$pop <- 1e6
  sim$count <- rpois(nrow(sim), exp(lr) * 10)
  std <- stats::setNames(rep(1, 8), ages)
  cl <- jp_cluster(sim, year, count, pop, age, std, kmax = 2,
                   ci_method = "parametric")
  expect_equal(sum(cl$selection$selected), 1)
  expect_equal(cl$selection$clusters[cl$selection$selected], 2)
  expect_equal(unique(cl$clusters$cluster[cl$clusters$age_group == "40-49"]),
               "40-49 to 70-79")
})

test_that("plots build without errors", {
  fit <- joinpoint(us_mortality, year, rate, by = "cause", kmax = 2,
                   control = jp_control(n_resample = 100, aapc_last = 10))
  plots <- list(
    autoplot(fit), autoplot(fit, color_by = "trend", log_y = TRUE,
                            joinpoint_ci = TRUE, aapc = TRUE, labels = "ci"),
    autoplot(fit, layout = "overlay", labels = "none"),
    plot_models(fit, 2), plot_apc(fit), plot_aapc(fit), plot_selection(fit),
    plot_residuals(fit), plot_residuals(fit, "qq"),
    plot_residuals(fit, "acf")
  )
  for (p in plots) expect_silent(ggplot2::ggplot_build(p))
})

test_that("results export to text files", {
  fit <- joinpoint(subset(us_mortality, cause == "Cancer"), year, rate,
                   kmax = 1, ci_method = "parametric")
  dir <- tempfile()
  on.exit(unlink(dir, recursive = TRUE))
  paths <- jp_export(fit, dir)
  expect_true(all(file.exists(paths)))
  expect_true(any(grepl("APCs", paths)))
})

test_that("age adjustment rejects incomplete or inconsistent age data", {
  d <- expand.grid(year = 1:2, age = c("a", "b"))
  d$count <- c(10, 20, 30, 40)
  d$pop <- 1000
  std <- c(a = 1, b = 3)
  expect_equal(nrow(age_adjust(d, year, count, pop, age, std)), 2)
  expect_error(age_adjust(d[-1, ], year, count, pop, age, std), "every time")
  expect_error(age_adjust(rbind(d, d[1, ]), year, count, pop, age, std),
               "Duplicated")
  expect_error(age_adjust(d, year, count, pop, age, c(a = -1, b = 3)),
               "invalid")
  d$std <- c(1, 2, 3, 3)
  expect_error(age_adjust(d, year, count, pop, age, std), "varies")
})

test_that("identical cohorts are not declared different", {
  d <- data.frame(x = rep(1:12, 2), g = rep(c("a", "b"), each = 12), y = 1)
  f <- joinpoint(d, x, y, by = "g", kmax = 0, log = FALSE,
                 ci_method = "parametric")
  cmp <- jp_compare(f, n_perm = 99)
  expect_equal(cmp$p_value, 1)
  expect_equal(cmp$result, "Parallel")
  expect_true(all(is.na(model_selection(f)$criteria$r_squared)))
})

test_that("ambiguous age labels are refused in clustering", {
  d <- expand.grid(year = 2000:2003, age = c("<5", "5-9", "10-14"),
                   stringsAsFactors = FALSE)
  d$count <- 10
  d$pop <- 1000
  expect_error(jp_cluster(d, year, count, pop, age,
                          standard = c("<5" = 1, "5-9" = 1, "10-14" = 1),
                          kmax = 0), "ambiguous")
})
