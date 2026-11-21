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

test_that("a dominant hinge leaves the later joinpoints and profile exact", {
  x <- 1:30
  y <- 1e8 * pmax(x - 10, 0) + 10 * pmax(x - 20, 0) + .01 * sin(x)
  m <- one(L(data.frame(x, y), kmin = 2, kmax = 2, method = "bic"))$models[["2"]]
  cf <- utils::combn(3:28, 2)
  cf <- cf[, cf[2, ] - cf[1, ] >= 3, drop = FALSE]
  oracle <- apply(cf, 2, function(t) {
    sum(lm.fit(cbind(1, x, pmax(x - t[1], 0), pmax(x - t[2], 0)), y)$residuals^2)
  })
  expect_equal(m$tau, cf[, which.min(oracle)])
  expect_lt(abs(m$sse - min(oracle)), 1e-5)
  expect_lt(abs(min(m$profile) - min(oracle)), 1e-5)
})

test_that("a large covariate component does not change the selection", {
  x <- 1:30
  z <- sin(.7 * x)
  base <- 10 * pmax(x - 15, 0) + .1 * cos(x)
  for (w in c("none", "se", "full")) {
    fit <- function(a) {
      d <- data.frame(x, y = base + a * z, z, s = seq(.5, 1.5, length.out = 30))
      args <- list(d, quote(x), quote(y), log = FALSE, type = "other",
                   kmax = 1, method = "bic", ci_method = "parametric",
                   covariates = ~ z)
      if (w == "se") args$se <- quote(s)
      if (w == "full") args$vcov <- outer(d$s, d$s) * .5^abs(outer(x, x, "-"))
      one(do.call(joinpoint, args))
    }
    f0 <- fit(0)
    f1 <- fit(1e7)
    expect_equal(f1$k_selected, f0$k_selected, info = w)
    expect_equal(diff(f1$criteria$bic), diff(f0$criteria$bic),
                 tolerance = 1e-5, info = w)
  }
})

test_that("permutation evidence ignores a large covariate component", {
  skip_on_cran()
  x <- 1:20
  z <- sin(.7 * x)
  base <- 10 * pmax(x - 10, 0) + .1 * cos(x)
  fit <- function(a) {
    one(L(data.frame(x, y = base + a * z, z), kmax = 1,
          method = "permutation", covariates = ~ z,
          control = jp_control(n_perm = 199)))
  }
  f0 <- fit(0)
  f1 <- fit(1e7)
  expect_equal(f1$k_selected, f0$k_selected)
  expect_equal(f1$permutation$table$p_value, f0$permutation$table$p_value)
})

test_that("the level of the response does not decide whether a slope is zero", {
  x <- 1:20
  y <- 1e-6 * x + 1e-8 * sin(x)
  fit <- function(a) {
    apc(joinpoint(data.frame(x, y = y + a), x, y, log = FALSE, type = "other",
                  kmax = 0, method = "bic",
                  control = jp_control(n_resample = 499)))
  }
  a0 <- fit(0)
  a1 <- fit(1e8)
  expect_gt(a1$lower, 0)
  expect_equal(a1$p_value, a0$p_value)
  expect_true(a1$significant)
})

test_that("a joinpoint on an irregular grid is found", {
  x <- c(1:5, 1e7 * (1:4))
  y <- pmin(x, 3)
  f <- L(data.frame(x, y), kmin = 1, kmax = 1, method = "bic")
  expect_equal(joinpoints(f)$estimate, 3)
})

test_that("criteria and the search keep small SSEs beside large changes", {
  x <- 1:30
  f <- L(data.frame(x, y = 1e5 * pmax(x - 10, 0) + pmax(x - 21, 0) +
                      .01 * sin(x)), kmax = 2, method = "bic")
  bic <- model_selection(f)$criteria$bic
  expect_equal(bic[2:3], c(1.033161, -9.240771), tolerance = 1e-5)
  expect_equal(f$fits[[1]]$k_selected, 2)
  f <- L(data.frame(x, y = 1e7 * pmax(x - 10, 0) + pmax(x - 21, 0) +
                      .01 * sin(x)), kmin = 2, kmax = 2, method = "bic")
  expect_equal(joinpoints(f)$estimate, c(10, 21))
  expect_lt(model_selection(f)$criteria$sse, .01)
})

test_that("fitted values keep their precision far from the origin", {
  x <- 2 * (1:8)
  f <- L(data.frame(x, y = 4 + .4 * x), kmax = 0, shift = 1e16)
  expect_lt(max(abs(residuals(f))), 1e-8)
  x <- 1:8
  d <- data.frame(x, y = 7 + 3.7 * x + .05 * sin(x))
  expect_equal(fitted(L(d, kmax = 0, shift = 1e15)), fitted(L(d, kmax = 0)),
               tolerance = 1e-10)
})

test_that("a constant series has an exact zero trend", {
  a <- apc(L(data.frame(x = 1:7, y = 10, s = 2:8), se = s, kmax = 0))
  expect_equal(a$estimate, 0)
  expect_equal(a$p_value, 1)
  expect_false(a$significant)
  expect_true(a$lower <= 0 && a$upper >= 0)
})

test_that("an unidentified joinpoint gets the whole admissible interval", {
  j <- joinpoints(L(data.frame(x = 1:20, y = 3 + .7 * (1:20)), kmin = 1,
                    kmax = 1, method = "bic"))
  expect_equal(c(j$lower, j$upper), c(3, 18))
})

test_that("an aliased parameter gets no standard error", {
  x <- 1:15
  z <- pmax(x - 7, 0) + (x == 7)
  f <- L(data.frame(x, y = 10 + x + 5 * pmax(x - 7, 0) + .02 * sin(x), z),
         kmin = 1, kmax = 1, covariates = ~ z)
  me <- model_estimates(f)
  expect_true(is.na(me$se[me$parameter == "Slope 2"]))
  expect_true(is.na(me$se[me$parameter == "Slope 2 - Slope 1"]))
  expect_true(is.finite(me$se[me$parameter == "Slope 1"]))
})

test_that("likelihood candidates are compared at convergence", {
  d <- data.frame(x = 1:10, n = c(1000, 1, 1000, 100, 1000, 10, 1, 1, 1, 1),
                  p = 10^c(6, 6, 7, 7, 4, 3, 2, 3, 2, 1))
  m <- joinpoint(d, x, count = n, population = p, family = "poisson",
                 kmin = 1, kmax = 1, method = "bic")
  expect_equal(joinpoints(m)$estimate, 4)
  expect_equal(model_selection(m)$criteria$deviance, 10988.34,
               tolerance = 1e-6)
  # A nearly collinear covariate: the covariance comes from the QR factor.
  x <- 1:16
  d <- data.frame(x, z = x + 1e-7 * sin(x), p = 1000,
                  n = c(5, 8, 9, 9, 8, 11, 13, 8, 14, 11, 9, 16, 18, 13, 15, 14))
  f <- joinpoint(d, x, count = n, population = p, family = "poisson",
                 kmax = 0, covariates = ~ z)
  ref <- summary(stats::glm(n ~ x + z, stats::poisson(), data = d,
                            offset = log(p / 1e5)))$coefficients[, 2]
  expect_equal(unname(sqrt(diag(vcov(f)))[2:3]), unname(ref[2:3]),
               tolerance = 1e-4)
})
