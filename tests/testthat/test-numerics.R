L <- function(d, ...) {
  joinpoint(d, x, y, log = FALSE, type = "other", ci_method = "parametric",
            ...)
}
one <- function(f) f$fits[[1]]

test_that("the grid search finds the least squares joinpoint", {
  # A large common trend and an irregular time grid.
  d <- data.frame(x = 1:12)
  d$y <- 1e6 * d$x + .03 * pmax(d$x - 7, 0) + rep(c(1e-4, -1e-4), 6)
  expect_equal(joinpoints(L(d, kmin = 1, kmax = 1))$estimate, 7)
  x <- c(1:5, 1e5 * (1:4))
  f <- L(data.frame(x, y = pmin(x, 3)), kmin = 1, kmax = 1, method = "bic")
  expect_equal(joinpoints(f)$estimate, 3)
  x <- 1:20
  y <- 1e8 * x + 10 * pmax(x - 10, 0) + .1 * sin(x)
  m <- one(L(data.frame(x, y), kmin = 1, kmax = 1, method = "bic"))$models[[1]]
  oracle <- vapply(3:18, function(t) {
    sum(lm.fit(cbind(1, x - mean(x), pmax(x - t, 0)), y)$residuals^2)
  }, numeric(1))
  expect_equal(m$tau, (3:18)[which.min(oracle)])
  expect_equal(m$sse, min(oracle), tolerance = 1e-5)
})

test_that("a common linear trend does not change the selection", {
  x <- 1:20
  y <- 10 * pmax(x - 10, 0) + .1 * sin(x)
  fit <- function(r) one(L(data.frame(x, y = r), kmax = 1, method = "bic"))
  f0 <- fit(y)
  f1 <- fit(y + 1e6 * x)
  expect_equal(f0$k_selected, 1)
  expect_equal(f1$k_selected, 1)
  expect_equal(diff(f1$criteria$bic), diff(f0$criteria$bic), tolerance = 1e-5)
})

test_that("WBIC selection does not depend on the origin of x", {
  set.seed(3)
  x <- 1:30
  d <- data.frame(x, y = 5 + .5 * x + .15 * pmax(x - 15, 0) + rnorm(30, 0, .7))
  f0 <- one(L(d, kmax = 2, method = "wbic"))
  f1 <- one(L(d, kmax = 2, method = "wbic", shift = 1e10))
  expect_equal(f1$criteria$weight, f0$criteria$weight, tolerance = 1e-5)
  expect_equal(f1$criteria$wbic, f0$criteria$wbic, tolerance = 1e-5)
  expect_equal(f1$k_selected, f0$k_selected)
})

test_that("the modified BIC shifts with the units of y as defined", {
  set.seed(8)
  x <- 1:30
  y <- 5 + .5 * x + .3 * pmax(x - 15, 0) + rnorm(30, 0, .5)
  crit <- function(a) one(L(data.frame(x, y = a * y), kmax = 1,
                            method = "mbic"))$criteria$mbic
  # Rescaling y by a moves MBIC(1) - MBIC(0) by -(2 / n) log(a).
  expect_equal(diff(crit(1000)) - diff(crit(1)), -2 * log(1000) / 30,
               tolerance = 1e-8)
})

test_that("an infeasible larger model leaves the smaller ones", {
  d <- data.frame(x = 1:5, y = c(2, 2.3, 2.4, 3.8, 5), z = pmax(1:5 - 3, 0))
  f0 <- L(d, kmax = 0, method = "bic", covariates = ~ z)
  expect_warning(f1 <- L(d, kmax = 1, method = "bic", covariates = ~ z),
                 "too few observations")
  expect_equal(one(f1)$k_selected, 0)
  expect_equal(fitted(f1), fitted(f0))
  expect_error(suppressWarnings(L(d, kmin = 1, kmax = 1, method = "bic",
                                  covariates = ~ z)), "lower `kmin`")
  # Enough observations, but every hinge equals the covariate.
  d <- data.frame(x = 1:9, y = c(2, 2.3, 2.4, 3.8, 5, 6, 7.4, 8, 9.1))
  d$z <- pmax(d$x - 5, 0)
  expect_error(L(d, kmin = 1, kmax = 1, method = "bic", covariates = ~ z,
                 control = jp_control(min_obs_end = 4)), "rank deficient")
})

test_that("too few observations lower kmax but not kmin", {
  d <- data.frame(x = 1:5, y = c(1, 3, 2, 1, 2))
  expect_warning(f <- L(d, continuous = TRUE, kmax = 2),
                 "too few observations")
  expect_length(coef(f, type = "hinge"), 3L)
  expect_error(suppressWarnings(L(d, continuous = TRUE, kmin = 2, kmax = 2)),
               "lower `kmin`")
})

test_that("inference needs residual degrees of freedom", {
  f <- joinpoint(data.frame(x = 1:2, y = c(1, 3)), x, y, log = FALSE,
                 control = jp_control(n_resample = 50))
  expect_true(is.na(apc(f)$p_value))
  f <- L(data.frame(x = 1:8, y = c(3, 7, 2, 8, 1, 9, 2, 8)), kmax = 2,
         method = "bic",
         control = jp_control(min_obs_end = 1, min_obs_between = 0))
  expect_true(all(model_selection(f)$criteria$df > 0))
})

test_that("a one-observation segment has no intercept or slope SE", {
  x <- 1:12
  y <- 1 + x + .1 * sin(x) + 3 * pmax(x - 2, 0)
  f <- L(data.frame(x, y), kmin = 1, kmax = 1, method = "bic",
         control = jp_control(min_obs_end = 1))
  est <- model_estimates(f)
  expect_true(is.na(est$se[est$parameter == "Intercept 1"]))
})

test_that("Monte Carlo ties count whatever the sign of the statistic", {
  mc <- joinpoint:::jp_mc_pvalue
  expect_equal(mc(-1, function(n) rep(-1, n), 99)[["p"]], 1)
  expect_equal(mc(2, function(n) rep(2, n), 99)[["p"]], 1)
  expect_equal(mc(Inf, function(n) rep(5, n), 99)[["p"]], 0.01)
})

test_that("seeds near the integer limit work for every cohort", {
  expect_false(is.na(joinpoint:::jp_cohort_seed(1500000000L, 2L, "vary")))
})
