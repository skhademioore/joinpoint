d <- subset(us_mortality, cause == "Stroke")

test_that("selection methods run and pick a model within bounds", {
  for (m in c("wbic", "bic", "bic3", "mbic", "dds")) {
    fit <- joinpoint(d, year, rate, kmax = 3, method = m,
                     ci_method = "parametric")
    k <- fit$fits[[1]]$k_selected
    expect_true(k >= 0 && k <= 3)
  }
})

test_that("empirical quantile intervals contain the estimates", {
  fit <- joinpoint(d, year, rate, kmax = 2,
                   control = jp_control(n_resample = 200))
  a <- apc(fit)
  expect_true(all(a$lower <= a$estimate & a$estimate <= a$upper))
  expect_true(all(is.na(a$statistic)))
  j <- joinpoints(fit)
  expect_true(all(j$lower <= j$estimate & j$estimate <= j$upper))
})

test_that("MADWD lowers the number of joinpoints", {
  fit0 <- joinpoint(d, year, rate, kmax = 3, method = "bic",
                    ci_method = "parametric")
  fit1 <- joinpoint(d, year, rate, kmax = 3, method = "bic",
                    ci_method = "parametric",
                    control = jp_control(madwd = 3))
  expect_lt(fit1$fits[[1]]$k_selected, fit0$fits[[1]]$k_selected)
  res <- model_selection(fit1)$madwd$result
  expect_equal(res, c(rep("Fail", length(res) - 1L), "Pass"))
})

test_that("error structures are applied", {
  pi <- subset(us_mortality, cause == "Influenza and Pneumonia")
  pi$deaths <- round(pi$rate * pi$population / 1e5)
  fit <- joinpoint(pi, year, count = deaths, population = population,
                   variance = "poisson", kmax = 2, ci_method = "parametric")
  expect_equal(fit$settings$variance, "poisson")
  fit_ar <- joinpoint(d, year, rate, kmax = 2, autocorrelation = "estimate",
                      ci_method = "parametric")
  expect_true(abs(fit_ar$fits[[1]]$rho) < 1)
  fit_rho <- joinpoint(d, year, rate, kmax = 2, autocorrelation = 0.3,
                       ci_method = "parametric")
  expect_true(is.matrix(fit_rho$fits[[1]]$S))
  expect_error(joinpoint(d, year, rate, variance = "se"), "standard errors")
})

test_that("zero counts get 0.5 under the documented conditions", {
  dd <- data.frame(t = 1:12, n = c(0, 3, 5, 2, 0, 6, 7, 4, 9, 8, 11, 10))
  fit <- joinpoint(dd, t, n, type = "count", variance = "poisson", kmax = 1,
                   ci_method = "parametric")
  expect_equal(sum(fit$fits[[1]]$data$corrected), 2)
  expect_equal(min(fit$fits[[1]]$data$r), 0.5)
})

test_that("exclusions keep observed points but drop them from fitting", {
  fit <- joinpoint(d, year, rate, kmax = 2, exclude = c(1997, 1998),
                   ci_method = "parametric")
  dat <- as.data.frame(fit)
  expect_equal(sum(!dat$included), 2)
  expect_equal(nobs(fit), nrow(d) - 2)
  expect_equal(max(apc(fit)$end), 1996)
})

test_that("jump and comparability ratio models recover a coding change", {
  set.seed(2)
  x <- 1980:2010
  y <- exp(3 + 0.02 * (x - 1980) + log(0.85) * (x >= 1999) +
             rnorm(length(x), sd = 0.01))
  dd <- data.frame(x = x, y = y)
  fj <- joinpoint(dd, x, y, kmax = 1, jump = 1998, method = "bic",
                  ci_method = "parametric")
  jt <- jump_estimates(fj, model = 0)
  expect_equal(jt$estimate[2], 0.85, tolerance = 0.02)
  fc <- joinpoint(dd, x, y, kmax = 1, jump = 1998, comparability_ratio = 0.85,
                  method = "bic", ci_method = "parametric")
  expect_equal(apc(fc)$estimate[1], 2.02, tolerance = 0.05)
  p <- predict(fc)
  expect_lt(abs(log(p$fit[1] / y[1])), 0.05)
})

test_that("linear models report slopes and average annual changes", {
  fit <- joinpoint(d, year, rate, log = FALSE, kmax = 1, method = "bic",
                   ci_method = "parametric")
  a <- apc(fit)
  m <- fit$fits[[1]]$models[[as.character(fit$fits[[1]]$k_selected)]]
  expect_equal(a$estimate, m$slopes)
  expect_false(is.null(aapc(fit)))
  expect_true(any(grepl("Slope", names(jp_table(fit)))))
})

test_that("cohorts with unusable data are skipped with a warning", {
  dd <- us_mortality
  dd$rate[dd$cause == "Cancer"][3] <- NA
  expect_warning(
    fit <- joinpoint(dd, year, rate, by = "cause", kmax = 1,
                     ci_method = "parametric"),
    "Cancer skipped"
  )
  expect_null(fit$fits[[2]])
  expect_error(joinpoint(dd[dd$cause == "Cancer", ], year, rate),
               "Missing values")
})

test_that("control settings and model limits are validated", {
  expect_error(jp_control(perm_seed = NA), "perm_seed")
  expect_error(jp_control(aapc_full = "yes"), "aapc_full")
  expect_error(jp_control(dds_cutoffs = c(200, 10)), "dds_cutoffs")
  tb <- subset(us_mortality, cause == "Tuberculosis")
  expect_error(joinpoint(tb, year, rate, kmin = -1), "whole numbers")
  expect_error(joinpoint(tb, year, rate, kmax = 1.5), "whole numbers")
})
